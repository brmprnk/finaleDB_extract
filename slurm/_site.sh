# Loads the site settings (slurm/site.env) for every slurm wrapper. Not an sbatch script.
# Sourced from the repository root (the wrappers cd to the submission directory). Stops with a message when site.env
# is missing or incomplete. Template: slurm/site.env.example.
if [ ! -f slurm/site.env ]; then
	echo "slurm/site.env not found: copy slurm/site.env.example to slurm/site.env and fill in your cluster's values" >&2
	exit 1
fi
. slurm/site.env
for _v in PREPROCESSING_DIR FRAGMENT_DIR GENOME; do
	eval "_val=\${${_v}:-}"
	if [ -z "${_val}" ]; then
		echo "slurm/site.env does not set ${_v}" >&2
		exit 1
	fi
done
RUN="${RUN:-}"
