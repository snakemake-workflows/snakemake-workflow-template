The `profiles/` directory can contain any number of subdirectories, each containing a `config.yaml` file with a [workflow-specific profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles):

`profiles/<specific_profile_name>/config.yaml`

The profile `profiles/default/config.yaml` will automatically be used by snakemake whenever you don't provide a workflow-specific profile via `--workflow-profile`.
This means that any resources or other (command line) arguments specified there, will implicitly be used when running this workflow.
Thus, as a workflow developer, only put configurations there that you expect to work in most environments, but which the users might want to tweak.
And for rule-specific resource setting, preferably provide generally applicable settings right in the rule definition, if necessary via [dynamic resource](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#dynamic-resources) specification---users can always override those in a profile, if they need to.

For any more specific profiles, use separate and clearly named subdirectories.
For example use `profiles/slurm/config.yaml` for a slurm-specific profile, or even something like `profiles/slurm_uni_xyz/config.yaml` for a particular institutional slurm compute cluster.

It is also good practice to add clear documentation comments for each entry in a (workflow) profile.
This should explain the respective entry, indicate what kind of values can be used and why a particular value or setting were chosen.
To this end, it is often helpful to provide links to relevant documentation pages, either from snakemake, a snakemake plugin or a specific cluster environment.

In general, we welcome pull requests for 3rd-party workflows you are working with to include such a profile for your specific compute environment.
