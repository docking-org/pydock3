#!/bin/bash

# req:
# EXPORT_DEST
# TMPDIR
# ARRAY_JOB_DOCKING_CONFIGURATIONS
# INPUT_DIR
# SKIP_IF_ALREADY_DONE

#optional:
# EXPORT_MOL2
# SLEEP_SECONDS_AFTER_COPYING_OUTPUT


# set default for unset vars
if [[ -z $EXPORT_MOL2 ]]; then
	EXPORT_MOL2=true
fi
if [[ -z $SLEEP_SECONDS_AFTER_COPYING_OUTPUT ]]; then
	SLEEP_SECONDS_AFTER_COPYING_OUTPUT=0
fi

# get scheduler job / task IDs
if ( ! [ -z $SLURM_ARRAY_JOB_ID ] ) && ( ! [ -z $SLURM_ARRAY_TASK_ID ] ); then
	SCHEDULER_NAME="slurm"
	JOB_ID=$SLURM_ARRAY_JOB_ID
	TASK_ID=$SLURM_ARRAY_TASK_ID
elif ( ! [ -z $JOB_ID ] ) && ( ! [ -z $SGE_TASK_ID ] ); then
	SCHEDULER_NAME="sge"
	#JOB_ID=$JOB_ID # already set by SGE
	TASK_ID=$SGE_TASK_ID
else
	echo "Scheduler job ID & task ID not found!"
	exit 1
fi

function log {
	echo "[$(date +%X)]: $@"
}

# log information about this job
log host=$(hostname)
log user=$(whoami)
log SCHEDULER_NAME=$SCHEDULER_NAME
log JOB_ID=$JOB_ID
log TASK_ID=$TASK_ID
log EXPORT_DEST=$EXPORT_DEST
log TMPDIR=$TMPDIR
log ARRAY_JOB_DOCKING_CONFIGURATIONS=$ARRAY_JOB_DOCKING_CONFIGURATIONS
log INPUT_DIR=$INPUT_DIR
log EXPORT_MOL2=$EXPORT_MOL2
log SLEEP_SECONDS_AFTER_COPYING_OUTPUT=$SLEEP_SECONDS_AFTER_COPYING_OUTPUT

# validate required environmental variables
for var in EXPORT_DEST DOCKFILES TMPDIR ARRAY_JOB_DOCKING_CONFIGURATIONS INPUT_DIR; do
	if [[ -z var ]]; then
		echo "One or more of the following require environmental variables are not defined: EXPORT_DEST DOCKFILES TMPDIR ARRAY_JOB_DOCKING_CONFIGURATIONS INPUT_DIR"
		exit 1
	fi
done

# initialize all our important variables & directories
JOB_DIR=${TMPDIR}/$(whoami)/${SCHEDULER_NAME}_${JOB_ID}_${TASK_ID}
DOCKFILES_TEMP=$JOB_DIR/working/dockfiles

#
log JOB_DIR=$JOB_DIR
log DOCKFILES_TEMP=$DOCKFILES_TEMP

#
OUTPUT=${EXPORT_DEST}/${TASK_ID}
LOG_OUT=${TMPDIR}/${SCHEDULER_NAME}_${JOB_ID}_${TASK_ID}.out
LOG_ERR=${TMPDIR}/${SCHEDULER_NAME}_${JOB_ID}_${TASK_ID}.err

# create directories
mkdir -p $JOB_DIR/working
mkdir -p $DOCKFILES_TEMP

#
mkdir -p $OUTPUT
chmod -R 777 $OUTPUT

# copy dockfiles
awk "\$1==${TASK_ID}{for (j=2; j<=NF; j++) print \$j}" "$ARRAY_JOB_DOCKING_CONFIGURATIONS" | xargs -I {} cp {} "$DOCKFILES_TEMP"
echo "dockfiles: "
ls $DOCKFILES_TEMP

# get dock executable path from array job docking configurations file (dockexec is last column)
DOCKEXEC=$(awk -v var="$TASK_ID" '$1 == var {print $NF}' "$ARRAY_JOB_DOCKING_CONFIGURATIONS")
log DOCKEXEC=$DOCKEXEC

# validate that there is exactly one INDOCK file
num_indock_files=$(ls "$DOCKFILES_TEMP"/INDOCK* 2>/dev/null | wc -l)
if [ ! "$num_indock_files" -eq 1 ]; then
    echo "There must be exactly one INDOCK file! Either none or multilple found in $DOCKFILES_TEMP"
	  exit 1
else
    # change copied indock file's name to "INDOCK" for convenience
    for file in "$DOCKFILES_TEMP"/INDOCK*; do
        new_name="INDOCK"
        mv "$file" "$DOCKFILES_TEMP"/"$new_name"
        echo "Renamed file $file to $new_name"
    done
fi

#
mkdir $JOB_DIR/dockfiles
for f in $DOCKFILES_TEMP/*; do
	ln -s $f $JOB_DIR/dockfiles/$(basename $f)
done
rm $JOB_DIR/dockfiles/INDOCK

#
find $INPUT_DIR -name '*.db2*' | sort > $JOB_DIR/working/split_database_index

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
fix_indock $DOCKFILES_TEMP/INDOCK $JOB_DIR/dockfiles/INDOCK

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
# 1. move results/restart marker to $OUTPUT (if no restart marker, remove it from $OUTPUT if present)
# 2. move logs to $OUTPUT
# 3. remove the working directory
function cleanup {
	nout=$(ls $OUTPUT | grep OUTDOCK | wc -l)

	if [ $nout -ne 0 ] && ! [ -f $OUTPUT/restart ]; then
		log "Something seems wrong, my output is already full but has no restart marker. Removing items present in output and replacing with my results."
		rm $OUTPUT/*
		nout=0
		nlog=0
	fi

	cp -p $JOB_DIR/working/OUTDOCK $OUTPUT/OUTDOCK.$nout
	sleep 20  # necessary in order to prevent bug witnessed using DockOpt with Slurm on Gimel where OUTDOCK fails to appear by the time job has left queue

	if $EXPORT_MOL2; then
	  cp -p $JOB_DIR/working/test.mol2.gz $OUTPUT/test.mol2.gz.$nout
  fi
	cp -p $LOG_OUT $OUTPUT/$nout.out
	cp -p $LOG_ERR $OUTPUT/$nout.err

	if [ -f $JOB_DIR/working/restart ]; then
		mv $JOB_DIR/working/restart $OUTPUT/restart
	elif [ -f $OUTPUT/restart ]; then
		rm $OUTPUT/restart
	fi

	chmod -R 777 $OUTPUT  # TODO: is this necessary? try to remove

	rm -rf $JOB_DIR

  sleep $SLEEP_SECONDS_AFTER_COPYING_OUTPUT
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

