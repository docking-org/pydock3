#!/bin/bash

# req:
# EXPORT_DEST
# INPUT_SOURCE
# DOCKFILES
# DOCKEXEC
# SHRTCACHE
# LONGCACHE
# JOB_ID
# SGE_TASK_ID

function log {
	echo "[$(date +%X)]: $@"
}

# log information about this job
log host=$(hostname)
log user=$(whoami)
log EXPORT_DEST=$EXPORT_DEST
log INPUT_SOURCE=$INPUT_SOURCE
log DOCKFILES=$DOCKFILES
log DOCKEXEC=$DOCKEXEC
log SHRTCACHE=$SHRTCACHE
log LONGCACHE=$LONGCACHE
log JOB_ID=$JOB_ID
log SGE_TASK_ID=$SGE_TASK_ID
log df=$(df)

# validate required environmental variables
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC SHRTCACHE LONGCACHE JOB_ID SGE_TASK_ID; do
  if [[ -z var ]]; then
    echo "One or more of the following require environmental variables are not defined: EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC SHRTCACHE LONGCACHE JOB_ID SGE_TASK_ID"
    exit 1
  fi
done

# initialize all our important variables & directories
JOB_DIR=${SHRTCACHE}/$(whoami)/${JOB_ID}_${SGE_TASK_ID}
INPUT_TARBALL=$(sed "${SGE_TASK_ID}q;d" $EXPORT_DEST/joblist | awk '{print $1}')
EXTRACT_DIR=$JOB_DIR/input/$(basename $INPUT_TARBALL)
DOCKFILES_TEMP=$JOB_DIR/working/dockfiles

#
log JOB_DIR=$JOB_DIR
log INPUT_TARBALL=$INPUT_TARBALL
log EXTRACT_DIR=$EXTRACT_DIR
log DOCKFILES_TEMP=$DOCKFILES_TEMP

#
OUTPUT=${EXPORT_DEST}/$(sed "${SGE_TASK_ID}q;d" $EXPORT_DEST/joblist | awk '{print $2}')
LOG_OUT=${SHRTCACHE}/rundock_${JOB_ID}_${SGE_TASK_ID}.out
LOG_ERR=${SHRTCACHE}/rundock_${JOB_ID}_${SGE_TASK_ID}.err

# bring directories into existence
mkdir -p $JOB_DIR/working
mkdir -p $EXTRACT_DIR
mkdir -p $DOCKFILES_TEMP

#
mkdir -p $OUTPUT
chmod -R 777 $OUTPUT

# failsafe: delete old jobs that failed to delete themselves
#n_old_jobs=$(find $SHRTCACHE/$(whoami) -maxdepth 1 -type d -mmin +240 | wc -l) # job directories older than 4 hours old to be removed
#if [ $n_old_jobs -gt 0 ]; then
(
	flock -x 9
	find $SHRTCACHE/$(whoami) -maxdepth 1 -type d -mmin +360 | xargs rm -r 2>/dev/null
	flock -u 9
)9>/dev/shm/rundock_purge_$(whoami).lock
#	synchronize_all_but_first "purge_old" "find $SHRTCACHE/$(whoami) -maxdepth 1 -type d -mmin +240 | xargs rm -r"
#fi

#
cp -a $DOCKFILES/. $DOCKFILES_TEMP

#
log "starting input extract"
tar -C $EXTRACT_DIR -xzf $INPUT_TARBALL
log "ending input extract"
mkdir $JOB_DIR/dockfiles
for f in $DOCKFILES_TEMP/*; do
	ln -s $f $JOB_DIR/dockfiles/$(basename $f)
done
rm $JOB_DIR/dockfiles/INDOCK
find $EXTRACT_DIR -name '*.db2*' | sort > $JOB_DIR/working/split_database_index

# tells this script to ignore SIGUSR1 interrupts
#trap '' SIGUSR1

if [ -f $OUTPUT/restart ]; then
	cp $OUTPUT/restart $JOB_DIR/working/restart
fi

function fix_indock {
	in=$1
	out=$2
	if [ -f $out ]; then
		rm $out
	fi

	while read -r line; do
		isfileline=$(echo $line | grep _file | awk '{print $1}')

		if ! [ -z $isfileline ] && ! [ $isfileline = "ligand_atom_file" ] && ! [ $isfileline = "output_file_prefix" ]; then
			label=$(echo $line | awk '{print $1}')
			filepath=$(echo $line | awk '{print $2}')
			filename=$(basename $filepath)
			basepath=""
			while ! [[ "$(basename $filepath)" == dockfiles* ]] && ! [ $(basename $filepath) = . ]; do
				[ -z $basepath ] && basepath=$(basename $filepath) || basepath=$(basename $filepath)/$basepath
				filepath=$(dirname $filepath)
			done
			if [ $(basename $filepath) = . ]; then
				# if the file path does not include some variation of "dockfiles", then we leave it as is
				basepath=$(echo $line | awk '{print $2}')
			else
				basepath=../dockfiles/$basepath
			fi
			labelcharc=$(printf $label | wc -c)
			nspaces=$((30-labelcharc))
			spaces=$(printf %-${nspaces}s " ")
			echo -e "$label$spaces$basepath" >> $out
		else
			echo "$line" >> $out
		fi
			
	done < $in
}

# only need to fix the INDOCK file once- don't want jobs to go all nutty because multiple processes are trying to mess with the INDOCK file
#(
#	flock -x 9
fix_indock $DOCKFILES_TEMP/INDOCK $JOB_DIR/dockfiles/INDOCK
#        flock -u 9
#)9>/dev/shm/dockfiles.${dockfileshash}.lock

log "starting dock"
pushd $JOB_DIR/working > /dev/null 2>&1

$DOCKEXEC $JOB_DIR/dockfiles/INDOCK &
dockpid=$!

function notify_dock {
	echo "notifying dock!"
	kill -10 $dockpid
}

trap notify_dock SIGUSR1

wait $dockpid
sleep 5 # bash script seems to jump the gun and start cleanup prematurely when DOCK is interrupted. This is stupid but effective at preventing this

# don't feel like editing DOCK src to change the exit code generated on interrupt, instead grep OUTDOCK for the telltale message
sigusr1=`tail OUTDOCK | grep "interrupt signal detected since last ligand- initiating clean exit & save" | wc -l`

log "finished! cleaning up"

# cleanup will:
# 1. remove extracted tarfiles
# 2. move results/restart marker to $OUTPUT (if no restart marker, remove it from $OUTPUT if present)
# 3. move logs to $OUTPUT
# 4. remove the working directory
function cleanup {
	rm -r $EXTRACT_DIR

	nout=$(ls $OUTPUT | grep OUTDOCK | wc -l)

	if [ $nout -ne 0 ] && ! [ -f $OUTPUT/restart ]; then
		log "Something seems wrong, my output is already full but has no restart marker. Removing items present in output and replacing with my results."
		rm $OUTPUT/*
		nout=0
		nlog=0
	fi

	chmod -R 777 $OUTPUT

	cp -p $JOB_DIR/working/OUTDOCK $OUTPUT/OUTDOCK.$nout
	cp -p $JOB_DIR/working/test.mol2.gz $OUTPUT/test.mol2.gz.$nout
	cp -p $LOG_OUT $OUTPUT/$nout.out
	cp -p $LOG_ERR $OUTPUT/$nout.err

	if [ -f $JOB_DIR/working/restart ]; then
		mv $JOB_DIR/working/restart $OUTPUT/restart
	elif [ -f $OUTPUT/restart ]; then
		rm $OUTPUT/restart
	fi

	rm -rf $JOB_DIR
}

popd > /dev/null 2>&1

if [ $sigusr1 -ne 0 ]; then
	echo "s_rt limit reached!"
	[ -z $SKIP_CLEANUP ] && cleanup
	exit 0
else
	[ -z $SKIP_CLEANUP ] && cleanup
	exit 0
fi

