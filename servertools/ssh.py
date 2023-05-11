import paramiko
import re
import warnings
from mse.system.directory import Directory
import stat
import os

#TODO: A single class that can do what all the three classes are doing below
#TODO: Restructure the classes  in this,
#TODO: combine it with Jochen's module


class ShellHandler:

    def __init__(self, host=None, user=None, port=22):

        config = paramiko.config.SSHConfig.from_file(open(os.path.expanduser("~/.ssh/config")))
        hostdata = config.lookup(host)
        hostdata["username"] = hostdata.pop("user")
        if not hostdata:
            assert host and user
            hostdata = {"hostname": host, "username": user, "port": port}

        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

        self.ssh.connect(**hostdata)

        self.ssh.__hostname__ = hostdata["hostname"]
        self.ssh.__username__ = hostdata["username"]
        channel = self.ssh.invoke_shell()
        self.channel = channel
        self.stdin = channel.makefile('wb')
        self.stdout = channel.makefile('r')

    def __del__(self):
        if hasattr(self,"ssh"):
            self.close()

    def close(self):
        self.ssh.close()

    def execute(self, cmd):
        """

        :param cmd: the command to be executed on the remote computer
        :examples:	execute('ls')
                    execute('finger')
                    execute('cd folder_name')
        """
        cmd = cmd.strip('\n')
        self.stdin.write(cmd + '\n')
        finish = 'end of stdOUT buffer. finished with exit status'
        echo_cmd = 'echo {} $?'.format(finish)
        self.stdin.write(echo_cmd + '\n')
        shin = self.stdin
        self.stdin.flush()

        shout = []
        sherr = []
        exit_status = 0
        for line in self.stdout:
            if str(line).startswith(cmd) or str(line).startswith(echo_cmd):
                # up for now filled with shell junk from stdin
                shout = []
            elif str(line).startswith(finish):
                # our finish command ends with the exit status
                exit_status = int(str(line).rsplit(maxsplit=1)[1])
                if exit_status:
                    # stderr is combined with stdout.
                    # thus, swap sherr with shout in a case of failure.
                    sherr = shout
                    shout = []
                break
            else:
                # get rid of 'coloring and formatting' special characters
                shout.append(re.compile(r'(\x9B|\x1B\[)[0-?]*[ -/]*[@-~]').sub('', line).
                             replace('\b', '').replace('\r', ''))

        # first and last lines of shout/sherr contain a prompt
        if shout and echo_cmd in shout[-1]:
            shout.pop()
        if shout and cmd in shout[0]:
            shout.pop(0)
        if sherr and echo_cmd in sherr[-1]:
            sherr.pop()
        if sherr and cmd in sherr[0]:
            sherr.pop(0)

        return shin, shout, sherr

    def clean_execute(self, cmd):
        stdin, stdout, stderr = self.execute(cmd)
        nstdout = []
        for i in stdout:
            if i.startswith(self.ssh.__username__) or self.ssh.__username__+"@" in i:
                continue
            elif "end of stdOUT buffer" in i:
                continue
            nstdout.append(i)
        # print also stdout and stderr

        return stdin, nstdout, stderr


class SFTPHandler:

    def __init__(self,
                 server:paramiko.client.SSHClient,
                 wrkpath=None):

        self.server = server
        self.sftp = server.open_sftp()
        self.sftp.chdir(wrkpath)
        self.basepath = wrkpath

    def __del__(self):
        self.server.close()

    def get_file(self, remotefile, localpath, safe=False):

        if safe:
            warnings.simplefilter("error")
        dir = os.path.dirname(localpath)
        if not dir:
            dir = "."
        if localpath in os.listdir(dir):
            warnings.warn("RemoteFile already exists on local", UserWarning)

        self.sftp.get(remotefile, localpath)

    def get_dirs(self, source, target):

        """
        Retrieves the contents of source directory(including subdirectories) to the target directory
        :param source: path on the server. The contents of which will be downloaded.(Path should be relative to sftp path)
        :param target: path on the local, to which contents will be downloaded to from the server.
        :return: None
        """
        # test if the source is a directory

        for item in self.sftp.listdir(source):
            filesattr = self.sftp.lstat(item)
            if stat.S_ISDIR(filesattr.st_mode):
                try:  # we only have (sub-)directories
                    os.mkdir(item)
                except FileExistsError:
                    pass
                self.get_dirs(source+"/"+item, os.path.join(target, item))
            else: # we have file simply copy it
                self.get_file(source+"/"+item, os.path.join(target, item))

    def put_file(self, localpath, remotepath, safe=False):

        if safe:
            warnings.simplefilter("error")
        # rdir = os.path.dirname(remotepath)
        # if not rdir:
        #     rdir = "."

        if remotepath in self.sftp.listdir(".") and remotepath != "./":
            warnings.warn("Localfile already exists on the local; if safe is set to false, The file will be over-write"
                          ,UserWarning)
        self.sftp.put(localpath, remotepath)

    def put_dirs(self, localpath, remotepath, safe=False):
        """
        Copies contents including subdirectories of a folder to a remote path. The local and remote paths must be
        relative to  working directory on local machine and remote server, respectively.
        :param localpath: path to a local folder
        :param remotepath: path on the remote, if not found it will be created in the sftp current working directory
        :param safe: If True, In case of over-write transfer will not occur and an error will be raise.
        :return: None
        """
        if not os.path.isdir(localpath):
            self.put_file(localpath, remotepath+"/"+localpath, safe=safe)

        else:
            for item in os.listdir(localpath):
                litem = os.path.join(localpath, item)
                print("%s" %litem)
                if os.path.isdir(litem):
                    print("%s is a directory" %litem)
                    # we have subfolder
                    try:
                        self.sftp.mkdir(remotepath+"/"+item)
                    except OSError: # means directory already exists
                        pass
                    self.put_dirs(litem, remotepath+"/"+item, safe=safe)
                else: #we have files here
                    self.put_file(litem, remotepath+"/"+item, safe=safe)

    def create_folder(self, path:str):
        self.sftp.mkdir(path)

    def bash_rsync_retrieve(self, remotepath, localpath, dry=False):
        """
        wrapper around rsync utility to retrieve files.Take care of '/' at the end of first path(remotepath), which can change the
        final structure of the (sub)folders.
        :param remotepath:
        :param localpath:
        :return:
        """
        Directory.rsync(self.server.__hostname__+":%s" %remotepath, localpath, dry=dry)

    def bash_rsync_put(self, localpath, remotepath, dry=False):
        """
        wrapper around rsync untility to put files into a remotepath. The paths must be relative to working
        directory on local machine and remote server. Be careful of the backslah "/" at the end for localpath,
        It may lead nested directory of localpath image at the remotepath.
        :param localpath:
        :param remotepath:
        :param dry:
        :return:
        """

        Directory.rsync(localpath, self.server.__hostname__+":%s" % remotepath, dry=dry)


class Slurminterface:

    def __init__(self, server, wrkdirec="/work/scratch/"):

        self.server = server
        self.wrkdirec = wrkdirec
        self.status = []
        self.jobids = []
        self.path = []

    def submit_job(self, path, scriptname="submit"):

        self.server.execute("cd %s" %path) # dive into the directory

        _, stdout, stderr = self.server.clean_execute("sbatch %s" %scriptname)

        for line in stdout:

            if "submitted batch job" in line: # for slurm output message
                break

        jobid = line.rsplit(maxsplit=1)[1] # the last digit contains the jobid

        self.path.append(path)
        self.jobids.append(jobid)

        self.server.execute("cd %s" %self.wrkdirec)
        return jobid

    def allsubmit_jobs(self, paths, scriptname="submit"):

        for p in paths:
            self.submit_job(p, scriptname=scriptname)

    def cancel_job(self, jobid):

        _, stdout, _ = self.server.clean_execute("scancel %s" %jobid)

    def cancel_alljobs(self):

        for i in self.jobids:
            self.cancel_job(i)

        self.show_jobs()

    def show_jobs(self):

        _, status, _ =self.server.clean_execute("squeue")

    def get_status(self, jobid):

        stdin, stdout, stderr = self.server.clean_execute("squeue -j %s" %jobid)
        print("".join(stdout))

        if not stderr:
            print("".join(stderr))

        status = {}
        for line in stdout:
            if line.startswith("-----"):
                continue
            else:
                st = [s for i, s in enumerate(line.split()) if i in [0, 2]]
                status[st[0]] = st[-1]

        self.status.append(status)

        return status

    def rstatus(self):

        for i in self.jobids:
            self.get_status(i)


    # def retrieve_job(self):
    #
    # 	sftp = self.SFTPHandler(self.server, self.wrkdirec)
    #
    # 	if "COMPLETED" in self.state or "CANCELED" in self.state:
    # 		sftp.retrieve(self.path, "%s" %self.wrkdirec)
    #
    # 	else:
    # 		print("%s is still in state %s" %(self.jobid, self.state))

