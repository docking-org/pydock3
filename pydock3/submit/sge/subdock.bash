#!/bin/bash
# req:
# EXPORT_DEST
# INPUT_SOURCE
# DOCKFILES
# DOCKEXEC
# SHRTCACHE
# LONGCACHE

RUNDOCK_PATH=$(dirname $0)/rundock.bash

failed=0
function exists {
	env_name=$1
	desc=$2
	if [ -z "${!env_name}" ]; then
		echo "expected env arg: $env_name"
		echo "arg description: $desc" 
		failed=1
	fi
}

function exists_warning {
	env_name=$1
	desc=$2
	default=$3
	if [ -z "${!env_name}" ]; then
		echo "optional env arg missing: $env_name"
		echo "arg description: $desc"
		echo "defaulting to $default"
		export $env_name="$default"
	fi
}

exists EXPORT_DEST "nfs output destination for OUTDOCK and test.mol2.gz files"
exists INPUT_SOURCE "nfs directory containing one or more .db2.tgz files OR a file containing a list of db2.tgz files"
exists DOCKFILES "nfs directory containing dock related files and INDOCK configuration for docking run"
exists DOCKEXEC "nfs path to dock executable"
exists SHRTCACHE "temporary local storage for job files"
exists LONGCACHE "longer term storage for files shared between jobs"
exists_warning QSUB_ARGS "additional arguments to qsub job submit command" ""
exists_warning QSUB_EXEC "qsub executable" "qsub"
exists_warning SGE_SETTINGS "script of SGE settings to source before using qsub" ""

if [ $failed -eq 1 ]; then
	echo "exiting with error"
	exit 1
fi

ligand_atom_file=$(grep -w ligand_atom_file $DOCKFILES/INDOCK | awk '{print $2}')
if [ "$ligand_atom_file" != "split_database_index" ]; then
	echo "ligand_atom_file option in INDOCK is [$ligand_atom_file], should be [split_database_index]. Please change and try again."
	exit 1
fi

if [ -w $DOCKFILES ]; then
	cat $DOCKFILES/* | sha1sum | awk '{print $1}' > $DOCKFILES/.shasum
fi

mkdir -p $EXPORT_DEST
n=1
njobs=0
warned=
if [ -d $INPUT_SOURCE ]; then
	get_input_cmd="echo $(find $INPUT_SOURCE -name '*.db2.tgz' | sort)"
else
	get_input_cmd="cat $INPUT_SOURCE"
fi

for input in $($get_input_cmd); do
	if ! [ -f $EXPORT_DEST/$n/OUTDOCK.0 ]; then
		echo $input $n
		njobs=$((njobs+1))
	elif ! [ -f $EXPORT_DEST/$n/test.mol2.gz.0 ]; then
		rm $EXPORT_DEST/$n/*
		echo $input $n
		njobs=$((njobs+1))
	elif [ -f $EXPORT_DEST/$n/restart ]; then
		echo $input $n
		njobs=$((njobs+1))
	fi
	n=$((n+1))
done > $EXPORT_DEST/joblist
n=$((n-1))

echo "submitting $njobs out of $n input files. $((n-njobs)) already complete"

var_args=
for var in EXPORT_DEST INPUT_SOURCE DOCKFILES DOCKEXEC QSUB_EXEC SGE_SETTINGS SHRTCACHE LONGCACHE; do
	[ -z "$var_args" ] && var_args="-v $var=${!var}" || var_args="$var_args -v $var=${!var}"
done

echo $var_args

echo "sleeping for 5 seconds before submitting, ctrl-C to stop"
if ! [ -z $SGE_SETTINGS ]; then
	source $SGE_SETTINGS
fi

sleep 5
echo "submitting jobs"
$QSUB_EXEC $var_args $QSUB_ARGS -cwd -S /bin/bash -q !gpu.q -t 1-$njobs $RUNDOCK_PATH
