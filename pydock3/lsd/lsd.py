import subprocess as sp
import sys, os
import argparse
import oyaml as yaml

from abc import ABC, abstractmethod

from pydock3.files import DirBase, FileBase
from pydock3.files import S3File, File

def get_fileobj_from_path(path, dir=False):
    if path.startswith("s3://"):
        if not dir:
            return S3File(path)
        else:
            return S3Dir(path)
    else:
        return File(path)

class LSDBase:

    def __init__(self, basedir : DirBase, dockfiles : DirBase, sdi):

        self.basedir = basedir
        self.dockfiles = dockfiles

    def initialize(self):

        DirBase.copy(self.dockfiles, self.basedir)
        FileBase.copy(self.sdi, self.basedir.file('sdi'))

    def completion_status(self, checksdi : set[int], checktop : bool) -> dict[str, set], dict[str, bool]:

        results_dict = {
            "lsd": {"success":set(), "failed":set(), "inactive":set()},
            "top": {"success":False, "failed":False, "inactive":False}
        }

        for f in (self.basedir + '/lsd').listall():
            sdi_index = int(f.dir().basename())
            if not sdi_index in checksdi:
                continue

            if f.path.endswith('.lsd_status'):
                if f.size() == 1:
                    results_dict["top"]["success"].add(sdi_index)
                else:
                    results_dict["top"]["failed"].add(sdi_index)

        results_dict["lsd"]["inactive"].update(
            checksdi.difference(
                results_dict["lsd"]["success"].union(results_dict["lsd"]["failed"])
            )
        )

        top_status  = (self.basedir.file('/.poses_status'))

        if top_status.exists():
            if top_status.size() == 1:
                results_dict["top"]["success"] = True
            else:
                results_dict["top"]["failed"]  = True
        else:
            results_dict["top"]["inactive"] = True

        return results_dict["lsd"], results_dict["top"]

class LSDQueue:

    def __init__(self, queue : JobQueue, lsdbase : LSDBase):
        self.queue = queue


    def submission_status(self) -> dict[str, set], dict[str, bool]:
        jobs = self.queue.get_jobs(jobfilter=lambda job: job.name.startswith(self.lsdbase.get_id()))

        results_dict = { 
            "lsd": {"running":set(), "submitted":set()} 
            "top": {"running":False, "submitted":False}
        }

        # <id>-<type>-<addtl>
        for job in jobs:
            lsd_id, lsd_job_type, addtl = job.name.split('-')

            if lsd_job_type == "lsd":
                lsd_sdi_offset = int(addtl)
                job_sdi_index = job.task_id + lsd_sdi_offset
                if job.running:
                    results_dict["lsd"]["running"].add(job_sdi_index)
                else:
                    results_dict["lsd"]["submitted"].add(job_sdi_index)

            elif lsd_job_type == "top":
                if job.running:
                    results_dict["top"]["running"] = True
                else:
                    results_dict["top"]["submitted"] = True

            elif lsd_job_type == "lsd_cont":
                rem_sdi_offset = int(addtl)
                rem_sdi = range(rem_sdi_offset, self.lsdbase.sdi_length())

                results_dict["lsd"]["submitted"].update(rem_sdi)

        return results_dict["lsd"], results_dict["top"]




    

def status_lsd(lsdbase, lsdqueue):

    all_jobs_lsd = set(range(lsdbase.sdi_length()))

    # let the queue figure out what jobs are currently in-transit
    sub_status = lsdqueue.submission_status(all_jobs_lsd)

    not_submitted_lsd = all_jobs.difference( sub_status["lsd"]["running"].union(sub_status["lsd"]["submitted"]) )
    not_submitted_top = not (sub_status["top"]["running"] or sub_status["top"]["submitted"])

    # the lsdbase (responsible for files) will be responsible for determining if a job has failed or succeeded (or isn't active yet)
    # technically in some cases the queue is also privy to this information, for example on aws batch 
    # but in the case of aws tracking job result status, those records are eventually purged, thus checking the file system is more reliable
    com_status = lsdbase.completion_status(not_submitted_lsd, not_submitted_top)

    # placeholder
    print(sub_status, com_status)

if __name__ == "__main__":

    mode = sys.argv[1]
    cfg  = sys.argv[2]

    cfg = yaml.load(cfg, Loader=yaml.Loader)

    queue, base = parse_config(cfg)

    if mode == "init":
        init_lsd(cfg)
    
    elif mode == "run":
        run_lsd(cfg)

    elif mode == "status":
        status_lsd(cfg)