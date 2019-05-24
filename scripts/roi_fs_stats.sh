
startdir=$(pwd)
source ~/PD_motor/global_scripts/setUpFreesurfer.sh

cd $SUBJECTS_DIR

for sub in */ ; do
	if [ "${sub:0:1}" == "0" ]; then
		echo Label file: $sub/label/lh.sensmotor.label
		echo Output file:  ~/PD_motor/rest_ec/fs_stats/${sub:0:4}.lh.sensmotor.stats
		mris_anatomical_stats \
			-l $sub/label/lh.sensmotor.label \
			-f ~/PD_motor/rest_ec/fs_stats/${sub:0:4}.lh.sensmotor.stats \
			$sub \
			lh
	fi
	
done

cd $startdir
