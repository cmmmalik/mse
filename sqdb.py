import os
import datetime

from Simtools.Jobs.job import Corejob
from Database.sqlite.sqlite import Job_sqlite

filename = os.path.expanduser("~/.HPCtools/Sim.db")


def checkarg(method):

    def deco(sq:Simsqlite, **kwargs):
        for k in kwargs:
            if not k in sq.columns:
                show = ",".join(sq.columns)
                raise ValueError(f"Invalid argument {k}, allowed arguments are \n{show}")

        return method(sq, **kwargs)

    return deco


class Decorator:

    @staticmethod
    def initialdec(init):

        def initialize(job:Corejob, **kwargs):

            init(job, **kwargs)
            database = Simsqlite(filename=filename)
            database.add_entry(jobname=job.name,
                           crtime="{}".format(datetime.datetime.now()),
                           localdir=job.working_directory,
                           status=job.status)

        return initialize


class Simsqlite(Job_sqlite):

    @checkarg
    def add_entry(self, **kwargs):
        cursor = self.insert(tablename=self.table_name, **kwargs)
        return cursor.lastrowid

    @checkarg
    def update_entry(self, Id, **kwargs):
        Set = [f"{k}={v}" for k,v in kwargs.items()]
        Set = ",".join(Set)
        self.update(tablename=self.table_name, Set=Set, condition=f"id={Id}")

    def delete_entry(self, Id):
        self.delete(tablename=self.table_name,condition=f"id={Id}")
