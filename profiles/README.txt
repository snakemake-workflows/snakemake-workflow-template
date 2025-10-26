A 'profiles' directory might contain workflow resource configuration for a particular cluster or cloud instance.

We encourage to include a profile for your use case as:

`profiles/<your cluster or cloud instance>/config.yaml`

You may include a readme file next to the config.yaml file to point out pitfalls or other things to consider.

Feel free to open pull requests for 3rd party workflows you are working with to include such a profile! It might also be necessary to occasionally label certain rules of a particular workflow with the `localrules: <rule 1>, <rule 2>, ...` directive when workflow developers focussed on server execution during development, e.g. plotting and dowload rules.
