from pydock3.files import DirBase, FileBase

from pydock3.lsd.lsd import load_lsd

# checkpoint job-
# if there are poses that can be processed, submit a job to select top poses
# if there are too few lsd jobs queued, submit more (with corresponding checkpoints)
def cont_lsd(cfg, workdir, run, block, stride, priority)

    lsdbase, queue = load_lsd(cfg)
    queue.load_job_definition("lsd")
    queue.load_job_definition("cont_lsd")

    sdi = lsdbase.sdi(run, block)

    for orig_block in sdi.orig_blocks():
        blocklabel = lsdbase.block_label(block)
        blockdir   = lsdbase.dir("lsd").dir(blocklabel)
        blockdir.create()

    lsdjobname = f"lsd-r{run}-b{block}"
    cntjobname = f"cont_lsd-r{run}-b{block}"

    # implementing a ghetto file lock here
    # still vulnerable to a race condition (esp with >2 concurrent processes), but this narrows the vulnerable time frame to the second or two (max) it takes to write this file rather than the entire length of the job
    # why not use a real locking mechanism?
    # short answer: can't be bothered
    # long answer: the method to do this will vary depending on storage provider
    # for example s3 does not provide a simple mechanism for acquiring a lock on a file- 
    # given that this race condition being met is very unlikely, it's not worth the effort

    # the "proper" way to do this in the cloud would be to spin up some sort of database that can serve as a lock store and have them access that
    # so consider this a TODO
    lockfile = lsdbase.file("lsdcheckpt")
    timeout = 120
    while lockfile.exists():
        timeout -= 1
        time.sleep(1)
        if timeout <= 0:
            raise Exception("timed out waiting for checkpoint lock to clear!")
    # here! is where the race condition will be met
    with lockfile.open('w') as lf:
        lf.write('locked!')

    try:
        next_block, n_active_lsd_jobs = collect_lsd_job_statistics(run)
        if not next_block and n_active_lsd_jobs == 0:
            last_block = True

        while next_block and n_active_lsd_jobs < queue.max_queued:

            if last_running_lsd_block < lsdbase.sdi.num_blocks(run):

                sdi_next = lsdbase.sdi.get(run, next_block)

                lsdjob = queue.submit(
                    "lsd", array=len(sdi_next), priority=priority,
                    params={
                        "LSD_CFG_PATH":lsdbase.cfg_path,
                        "LSD_WORKDIR":workdir,
                        "LSD_RUN":run,
                        "LSD_BLOCK":block
                    }
                )

                queue.submit(
                    "cont_lsd", deps=lsdjob, priority=priority, 
                    params={
                        "LSD_CFG_PATH":lsdbase.cfg_path, 
                        "LSD_WORKDIR":workdir,
                        "LSD_RUN": run, 
                        "LSD_BLOCK":block+stride, 
                        "LSD_SUBMIT_STRIDE":stride, 
                        "LSD_SUBMIT_PRIO":priority
                    }
                )

                next_block, n_active_lsd_jobs = collect_lsd_job_statistics(run)

        free_poses, existing_top_job_id, existing_top_job_result = get_top_stats()
        total_free_poses = sum([fp.nposes for fp in free_poses])

        posesthreshold = lsdbase.top.n * 5
        if last_block or freeposes > posesthreshold:

            topjob = queue.submit(
                "top", deps=existing_top_job_id, 
                params={
                    "LSD_CFG_PATH":lsdbase.cfg_path,
                    "LSD_WORKDIR":workdir,
                    "LSD_TOP_POSES_SDI":newposesfile,
                }
            )

    except:
        traceback.print_exc()
        return 1
    # make sure to clean up this "lock" file!
    finally:
        lockfile.rm()
    return 0

if __name__ == "__main__":

    cfgpath     = os.getenv("LSD_CFG_PATH")
    workdir     = os.getenv("LSD_WORKDIR")
    run         = os.getenv("LSD_RUN")
    block       = os.getenv("LSD_BLOCK")
    stride      = os.getenv("LSD_SUBMIT_STRIDE")
    priority    = os.getenv("LSD_SUBMIT_PRIO")

    sys.exit(cont_lsd(cfgpath, run, block, stride, priority))