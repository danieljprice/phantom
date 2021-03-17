Long-term archiving of your phantom calculations
==================================================================
One of the biggest headaches we have with shared supercomputer projects
is that inevitably somebody fills whatever disk quota was allocated,
and the project halts for everyone. A true tragedy of the commons. To solve this, shift your data somewhere more permanent.

Data curation
-------------
For calculations with phantom that have been published in a paper,
best practice is to upload the **entire calculation including .in and
.setup files, .ev files and all dump files in a public repository**.

See for example a dataset from Mentiplay et al. (2020) using figshare: `https://doi.org/10.6084/m9.figshare.11595369.v1`_

Or this example from Wurster, Bate & Price (2018) in the University of Exeter repository: `https://doi.org/10.24378/exe.607`_

Archiving your data to Google Drive using rclone
------------------------------------------------
You can use rclone to copy data from a remote cluster or supercomputing facility to Google Drive. For universities with institutional subscriptions, this provides almost unlimited storage.

Set this up by logging into your supercomputer and typing::

   $ rclone config
   No remotes found - make a new one
   n) New remote
   s) Set configuration password
   q) Quit config

   name> dan-google-drive

   Storage> drive

   Google Application Client Id
   See https://rclone.org/drive/#making-your-own-client-id for how to create your own.

   client_id>   (enter client id you got from the instructions)
   client_secret>  (enter client secret you got from the instructions)

   scope> drive
   root_folder_id> (leave this blank)

   Edit advanced config? (y/n) n

   Use auto config?
   y/n> n

   Please go to the following link: https://accounts.google.com/o/oauth2/auth?access_type=offline&client_id=...
   (click on the link to approve)

   Configure this as a team drive?
   y/n> n

Check the above was successful by listing files on your remote drive using::

    $ rclone ls dan-google-drive:

To copy files to your google drive you can then use::

    $ rclone copy local_path remote_path

For example::

    $ rclone copy $HOME/runs/phantom/disc-test1 dan-google-drive:phantom/disc-test1

To sync/update an entire directory tree onto your google drive you can then use::

    $ rclone sync -i $HOME/runs dan-google-drive:runs

Other helpful information
--------------------------
- :doc:`General instructions for running on a remote cluster <running-clusters>`
- `rclone userguide <https://rclone.org>`_
