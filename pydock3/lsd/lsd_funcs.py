
class SDI(ABC):
    def __init__(self, basedir):
        self.basedir = basedir
        self.blocksize = 2000
        self.sdi_cache = {}
        self.sdi_cache_fifo = []
        self.cache_limit = int(16e6/self.blocksize/100)

    def get_entry_orig(self, index):
        pass

    def _get_cached_sdi(self, subset, block):
        return self.sdi_cache.get(f"s{subset}b{block}")

    def _set_cached_sdi(self, subset, block, data):
        # pop oldest entry when we reach the limit
        if len(self.sdi_cache) + 1 > self.cache_limit:
            pop_k = self.sdi_cache_fifo.pop()
            self.sdi_cache.pop(pop_k)
        self.sdi_cache[f"s{subset}b{block}"] = data
        self.sdi_cache_fifo.append()

    def get_entry(self, subset, block, index) -> int, str, DirBase:
        sdi_lines = self._get_cached_sdi(subset, block)
        if not sdi_lines:
            blockfile = self.sdidir_base.dir('s{subset}.d').file('b{block}.sdi')
            if not blockfile.exists():
                raise Exception("tried to access non-existent sdi block!")
            with blockfile.open('r') as bf:
                sdi_lines = bf.readlines()
            self.set_cached_sdi(subset, block, sdi_lines)
        orig_idx, entry_val = sdi_lines[index].split('\t')
        return orig_idx, entry_val, self._get_entry_outdir(orig_idx)

    def list_entries(self, subset, block) -> list[str]:
        sdi_lines = self._get_cached_sdi(subset, block)
        return sdi_lines

    def list_blocks(self, subset):
        block_i = 0
        block_f = self.sdidir_base.dir('s{subset}').file('b{block_i}.sdi')
        block_l = []
        while block_f.exists():
            block_l.append(block_i)
            block_i += 1
            block_f = self.sdidir_base.dir('s{subset}').file('b{block_i}.sdi')
        return block_l

    def list_subsets(self):
        subset_i = 0
        subset_f = self.sdidir_base.dir('s{subset_i}').file('b0.sdi')
        subset_l = []
        while subset_f.exists():
            subset_l.append(subset_i)
            subset_i += 1
            subset_f = self.sdidir_base.dir('s{subset_i}').file('b0.sdi')
        return subset_l
            
    def add_subset(self, subset_indices : list[int]):
        subset = len(self.list_subsets())
        orig_file = self.sdidir_base.file('orig')
        new_lines = []
        all_lines = []

        with orig_file.open('r') as of:
            all_lines = of.readlines()

        for idx in subset_indices:
            new_lines.append(f"{idx}\t{all_lines[idx]}")

        block = 0
        i = 0
        currfile = self.sdidir_base.dir('s{subset}.d').file('b{block}.sdi')
        currhndl = currfile.open('w')
        while i < len(new_lines):
            currhndl.write(new_lines[i])
            i += 1
            if i % self.blocksize == 0:
                block += 1
                currhndl.close()
                currfile = self.sdidir_base.dir('s{subset}.d').file('b{block}.sdi')
                currhndl = currfile.open('w')
        currhndl.close()
        return subset

    def basedir(self):
        return self.sdidir_base

    def subsetdir(self, subset):
        return self.basedir().dir('s{subset}.d')

    def blockdir(self, subset, block):
        return self.subsetdir(subset).dir('b{block}.sdi')

class ComputeQueue(ABC):
    @abstractmethod
    def ensure_path_is_accessible(fileobj):
        pass
    @abstractmethod
    def is_file_local(fileobj):
        pass
    @abstractmethod
    def submit(self, job, deps=None, array=None, params={})
        pass

class QBlasterBase(ABC):
    # help, of course- the most important command
    @abstractmethod
    def help(self):
        pass

    # top level submit and status for program
    @abstractmethod
    def submit(self, sdi_id=None):
        pass
    @abstractmethod
    def status(self, sdi_id=None):
        pass

    # add, subtract, and list standard input files for program
    @abstractmethod
    def add_sdi(self, sdi_path):
        pass
    @abstractmethod
    def rm_sdi(self, sdi_id) # sdi_id can be a checksum or real path to the file
        pass
    @abstractmethod
    def list_sdi(self):
        pass

    # get and set the compute queue- include checks to make sure files in base are accessible from queue
    @abstractmethod
    def set_queue(self, queue_name):
        pass
    @abstractmethod
    def get_queue(self, queue_name):
        pass

    # create an archive from the current base, deactivating it
    @abstractmethod
    def create_archive(self):
        pass
    # list all archives and their properties
    @abstractmethod
    def list_archives(self):
        pass
    # destroy an archive
    @abstractmethod
    def destroy_archive(self, archive_path):
        pass
    # destroy all files in this base
    @abstractmethod
    def destroy_base(self):
        pass

    @staticmethod
    def cp(self, p1, p2):
        f1 = fileobj_from_path(p1)
        f2 = fileobj_from_path(p2)
        FileBase.copy(p1, p2)

    @staticmethod
    def rm(self, p):
        f = fileobj_from_path(p)
        f.rm()

class JobStates(Enum):
    
    SUCCESS=0
    FAILURE=1
    PARTIAL=2
    INACTIVE=3
    SUBMITTED=4
    RUNNING=5
    UNKNOWN=6
    def from_string(s):
        s_dict = {
            "INACTIVE":JobStates.INACTIVE,
            "SUBMITTED":JobStates.SUBMITTED,
            "RUNNING":JobStates.RUNNING,
            "SUCCESS":JobStates.SUCCESS,
            "FAILURE":JobStates.FAILURE,
            "PARTIAL":JobStates.PARTIAL
        }
        if not s_dict.get(s):
            return JobStates.UNKNOWN
        else:
            return s_dict[s]

class LSDBase(QBlasterBase):
    def __init__(self, cfg, filebase : ParallelJobFileBase, queue : JobQueue):
        self.filebase = filebase
        self.sdibase = filebase.dir('sdi')
        self.jtype = "lsd"

    def _get_queue_stats(self, sdifilter=None, subsetfilter=None):
        sdim = {}
        qstats = {}

        for job in self.queue.listjobs():
            jtype, sdi_id, subset, block = j.name.split('-')
            subset = int(subset)
            block = int(block)

            if jtype != self.jtype:
                continue

            if sdifilter and sdi_id != sdifilter:
                continue
            if subsetfilter and subset != subsetfilter:
                continue

            sdi = sdim.get(sdi_id)
            if not sdi:
                sdi = self.get_sdi(sdi_id)
                sdim[sdi_id] = sdi
                qstats[sdi_id] = {}

            orig_idx = sdi.get_orig_idx(subset, block, j.task_id)
            qstats[sdi_id][orig_idx] = j.state

        return qstats

    def _get_sdi_stats(self, sdi, qstats=None):
        stats = {}
        if not qstats:
            qstats = self._get_queue_stats(sdifilter=sdi.id)
            qstats = qstats.get(sdi.id) or {}

        started_subsets = self._get_subsets_started(sdi)
        complete_subsets = self._get_subsets_complete(sdi)

        for subset in sdi.list_subsets():
            if subset in complete_subsets:

    def _get_subset_stats(self, sdi, subset, qstats=None):
        if not qstats:
            qstats = self._get_queue_stats(sdi, sdifilter=sdi.id, subsetfilter=subset)
            qstats = qstats.get(sdi.id) or {}

        ss_complete = True
        stats = self._load_subset_stats(sdi, subset)
        started_blocks = self._get_blocks_started(sdi, subset)
        complete_blocks = self._get_blocks_complete(sdi, subset)

        for block in sdi.list_blocks(subset):

            if (block in complete_blocks):
                continue

            if not (block in started_blocks):
                for orig_idx, entry in sdi.list_entries(subset, block):
                    qstats[orig_idx] = JobStates.SUBMITTED

            bk_complete = True
            for orig_idx, entry in sdi.list_entries(subset, block):
                if qstats.get(orig_idx):
                    bk_complete = False
                    continue
                elif stats.get(orig_idx):
                    continue
                else:
                    stats.update(self._get_entry_stats(sdi, subset, orig_idx))

            if not bk_complete:
                ss_complete = False
            else:
                self._mark_block_complete(sdi, subset, block)

        if ss_complete:
            self._mark_subset_complete(sdi, subset)

        self._save_subset_stats(sdi, subset, stats)

        return stats

    def _get_entry_stats(self, sdi, subset, orig_index):
        outdir = sdi.get_entry_outdir(index)
        statsf = outdir.file(f'{subset}-stats.json')
        if not statsf.exists():
            return JobStates.INACTIVE, {}
        with statsf.open('r') as statsf_obj:
            stats = json.loads(statsf_obj.read())
        # stats["jobstate"] should be one of SUCCESS, FAILURE, or PARTIAL
        state = JobStates.from_string(stats["jobstate"])
        del stats["jobstate"]
        return state, stats

    def _load_subset_stats(self, sdi, subset):
        statfile = sdi.subsetdir(subset).file('stats.json')
        return self._load_stats(statfile)

    def _save_subset_stats(self, sdi, subset, stats):
        statfile = sdi.subsetdir(subset).file('stats.json')
        return self._save_stats(statfile, stats)

    def _load_stats(self, statfile):
        if not statfile.exists():
            return {}
        with statfile.open('r') as stat_fo:
            stats = json.loads(stat_fo.read())
            return stats

    def _save_stats(self, statfile, stats):
        with statfile.open('w') as stats_fo:
            json.dump(stats, stats_fo)

    # functions for setting various completion markers
    # all of these assume a lock has been acquired for exclusive access to relevant files
    def _get_subsets_started(self, sdi):
        return self._load_completion(sdi.basedir().file('started'))

    def _get_subsets_complete(self, sdi):
        return self._load_completion(sdi.basedir().file('complete'))

    def _get_blocks_started(self, sdi, subset):
        return self._load_completion(sdi.subsetdir(subset).file('started'))

    def _get_blocks_complete(self, sdi, subset):
        return self._load_completion(sdi.subsetdir(subset).file('complete'))

    def _mark_block_started(self, sdi, subset, block):
        self._save_completion(sdi.subsetdir(subset).file('started'), block)

    def _mark_block_complete(self, sdi, subset, block):
        self._save_completion(sdi.subsetdir(subset).file('complete'), block)

    def _mark_subset_started(self, sdi, subset):
        self._save_completion(sdi.basedir().file('started'), subset)

    def _mark_subset_complete(self, sdi, subset):
        self._save_completion(sdi.basedir().file('complete'), subset)

    def _load_completion(self, cfile):
        with cfile.open('r') as cfile_o:
            data = cfile_o.read().split('\n')
            return [int(d) for d in data]

    def _save_completion(self, cfile, citem):
        with cfile.open('a') as cfile_o:
            cfile_o.write(f'{citem}\n')

    ### internal stats logic, used by jobs & status command
    # returns a string representing overall state, followed by a list of all jobs & their individual states
    # overall state can be:
    # state name  | allowed to submit?
    # inactive    | yes
    # running     | no
    # submitted   | no
    # failed      | yes
    # succeeded   | optional

    ### queue logic
    def set_queue(self, queue_name):
        queue_obj = load_queue(queue_name)
        if not queue_obj:
            print("Could not find queue!")
            return
        else:
            queue_obj.ensure_path_is_accessible(self.filebase)
            self.queue = queue
            
    def get_queue(self):
        

    ### dockfiles logic
    def set_dockfiles(self, dockfiles_path, force=False):
        curr = self.get_dockfiles()
        if force or not curr:
            df_src = dirobj_from_path(dockfiles_path)
            df_dst = self.filebase.dir("dockfiles")
            DirBase.copy(df_src, df_dst)
            print("Successfully copied over dockfiles")
        else:
            rough_status = self.get_rough_lsd_status()
            if rough_status == "submitted":
                print("Error! Active jobs are currently using dockfiles. Wait for them to finish or cancel them.")
            elif rough_status in ["partial", "complete"]:
                print("Warning! You are trying to change your dockfiles but some poses have already been produced from them.")
                print("Choose one of the following actions:")
                print("1. Archive existing work, initialize a fresh run with new dockfiles")
                print("2. Destroy existing work, initialize a fresh run with new dockfiles")
                print("3. Continue with existing run, just swap out dockfiles.")
                choice = input("[choose from 1, 2, or 3]: ")

                if choice == 1:
                    self.create_archive()
                    self.set_dockfiles(dockfiles_path)
                elif choice == 2:
                    self.destroy_base()
                    self.set_dockfiles(dockfiles_path)
                elif choice == 3:
                    self.set_dockfiles(dockfiles_path, force=True)

    def get_dockfiles(self):
        df = self.filebase.dir("dockfiles").file("INDOCK")
        if df.exists():
            return df.dir()
        else:
            return None

    ### archive/delete logic
    def list_archives(self):
        pass
    def create_archive(self):
        pass
    def destroy_base(self):
        pass

    ### SDI logic
    def add_sdi(self, sdi_path):

        new_sdi = SDI.from_path(sdi_file)
        if self.get_sdi(new_sdi.signature):
            print("sdi already exists!")
            return False
        else:
            new_sdi.source_path = sdi_path
            sdi_dest = self.filebase.dir("sdi").file(new_sdi.signature)
            new_sdi.save(sdi_dest)
        return new_sdi.signature

    def get_sdi(self, sdi_id):
        
        sdi_file = self.sdibase.file(sdi_id)
        if sdi_file.exists():
            sdi = SDI.from_file(sdi_file)
            return sdi
        else:
            return None

    def rm_sdi(self, sdi_id):

        sdi_file = self.sdibase.file(sdi_id)
        if sdi_file.exists():
            sdi = SDI.from_file(sdi_file)
            lsd_status, lsd_stats = get_lsd_stats(sdi.signature)
            if lsd_status != "inactive":
                print("")
            sdi.destroy()

    def list_sdi(self):
        pass

    ### submit logic
    def __submit_sdi(self, sdi_id, dockfilesname, submitlist=None):

        self.queue.acquire_lock(timeout=60, name=f"submit-lsd-{sdi_id}")

        try:
            sdi_blocks_status, sdi_flattened_status = self.get_lsd_stats(sdi_id)

            # you are not permitted (period) to submit a job that is already submitted/running
            banned_resubmit_states = ["submitted", "running"]
            resub_stats = {
                "inactive" : 0,
                "failed"   : 0,
                "partial"  : 0, # "partial" is an optional state for jobs like lsd that can partially complete
                "succeeded": 0
            }
            if not submitlist:
                default_resubmit_states = ["inactive", "failed"]
                sdi_ready_list = filter(lambda x: (x[1] in default_resubmit_states) and (x[1] not in banned_resubmit_states), sdi_flattened_status)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:x[1] == resub_type, sdi_ready_list))
                sdi_ready_list = [x[0] for x in sdi_ready_list]
                
            # if desired, you *may* resubmit a job that has already succeeded or is in some state other than submitted/running
            elif submitlist:
                sdi_ready_list = filter(lambda x: (flattened_status[x] not in banned_resubmit_states), submitlist)
                for resub_type in resub_stats:
                    resub_stats[resub_type] = len(filter(lambda x:flattened_status[x]==resub_type, submitlist))

            sdi = self.get_sdi(sdi_id)
            subset = sdi.add_subset(subset_indices=sdi_ready_list)

            # no actual submission happens yet- we merely mark it as such
            # this way our continuation submission has a record of all jobs to submit blocks for
            self.__mark_subset_submitted(sdi_id, subset)

        finally:

            self.queue.release_lock(name=f"submit-{sdi_id}")

    def submit(self, sdiname=None, confirm=False):

        if sdiname is None:
            for sdiname in self.listsdi():
                self.__submit_sdi(sdiname)
        else:
            self.__submit_sdi(sdiname)

        self.continue_submission()

    # continuation function- very important, does the actual submission to queue. "submit" just gathers statistics & creates new subsets/submit markers for continuation to use
    # continue does not operate on a single sdi/subset, it operates on the submission queue as a whole
    # ergo all sdi/subsets that have been submitted are *actually* submitted through this function
    # this allows all jobs to share the same queue, while not requiring a dedicated process that constantly exists and monitors for free space in the queue per sdi
    def continue_submission(self):

        self.queue.acquire_lock(timeout=60, name=f"continue")

        # psuedocode:
        # for sdi_id, subset, block in get_available_blocks_to_submit():
        #   queue.submit('lsd', array=blocksize(block), params={sdi_id, subset, block})
        #   mark_block_submitted(sdi_id, subset, block)

    ### status logic
    def __status_sdi(self, sdi):
        stats = []
        for run in sdi.runs():
            if run.stats:
                stats += run.stats
                continue
            for block in run.blocks():
                if block.stats:
                    stats += block.stats
                    continue
                for entry in block:
                    stats += entry.stats
    
    def status(self, sdiname=None):

        pass

    ###################
    # example session #
    ###################

    # qblaster ligands --base /nfs/home/mystuff

    # qblaster ifp

    # qblaster ecfp

    # qblaster lsd --base s3://mybucket/mylsdbase
    # > add-sdi z22://h17p200-h17p400
    # > add-sdi z22://h18p200-h18p400
    # > set-dockfiles ~/mydockfiles
    # > get-dockfiles
    #   +------------------------------------+----------+
    #   | path                               | checksum |
    #   +------------------------------------+----------+
    #   | s3://mybucket/mylsdbase/dockfiles  | 123fda   |
    #   +------------------------------------+----------+
    # > set-queue aws_ue1 # where are access keys set?
    # > submit
    # > set-dockfiles ~/otherdockfiles
    #   Error! Active jobs are currently using dockfiles. Wait for them to finish or cancel them.
    #
    # > set-queue slurm
    #   Error! Active jobs are currently using queue. Wait for them to finish or cancel them.
    #
    # > cancel
    #   Cancelling Jobs... Done!
    #
    # > set-dockfiles ~/otherdockfiles
    #   Warning! You are trying to change your dockfiles but some poses have already been produced from them.
    #   Choose one of the following actions:
    #   1. Archive existing work, initialize a fresh run with new dockfiles
    #   2. Destroy existing work, initialize a fresh run with new dockfiles
    #   3. Continue with existing run, just swap out dockfiles.
    #   [choose from 1, 2, or 3]: 1
    #   Archive just top poses? [y/N]: y
    #
    #   Creating archive @ s3://mybucket/mylsdbase/123fda.archive.zip
    #
    # > set-queue slurm
    #   Warning: base in "S3" is accessible but not local to queue "slurm", expect additional data transfer fees
    # 
    # > list-sdi
    #   +-----------------------+-----------+
    #   | path                  | checksum  |
    #   +-----------------------+-----------+
    #   | z22://h17p200-h17p400 | a342ve    |
    #   | z22://h18p200-h18p400 | 1usq3r    |
    #   +-----------------------+-----------+
    # > submit z22://h17p200-h17p400
    # > status
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | sdi                   | inactive  | submitted | running   | success   | failure   | poses     | total time  | poses data size | total data processed |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    #   | z22://h17p200-h17p400 | 0         | 9375      | 604       | 20        | 1         |
    #   | z22://h18p200-h18p400 | 10000     | 0         | 0         | 0         | 0         |
    #   +-----------------------+-----------+-----------+-----------+-----------+-----------+
    # > status --top
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | sdi                   | harvested poses | free poses      | top poses | poses path                                 | scores path                               |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    #   | z22://h17p200-h17p400 | 9503136         | 100324          | 300000    | s3://mybucket/mylsdbase/top/a342ve.mol2.gz | s3://mybucket/mylsdbase/top/a342ve.scores |
    #   | z22://h18p200-h18p400 | 0               | 0               | 0         | N/A                                        | N/A                                       |
    #   +-----------------------+-----------------+-----------------+-----------+--------------------------------------------+-------------------------------------------+
    # > cp 
    # > rm
    # > 