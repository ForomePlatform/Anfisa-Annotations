# Annotation Database Service

## Overview

For an efficient annotation process we need to store annotation sources locally
for the actual process. To annotate a single genetic variant, in most cases, we
can access all the information using a well-defined key, based on the genetic
coordinates of the variant and the DNA change. Hence, a key-value store seems to
be a perfect place to store the majority of the information. We use RocksDB as
the key-value store for annotations.  
                                
## Configuration Notes (for Linux)

> If the annotation service is being set up on a Linux host, 
> it is important to ensure that the number of files that can be 
> simultaneously opened is no less than 10,000. This can be checked 
> with the following Linux command:

    ulimit -S

> To change the value if it is too low, add a line to 
> `/etc/security/limits.conf`:

    <user>       	soft	nofile      	20000
> Additionally, in files `/etc/pam.d/common-session-noninteractive`
> and `/etc/pam.d/common-session`, add a line:

    session required pam_limits.so
> Finally, ensure, that in /etc/ssh/sshd_config there is a line:
    
    UsePAM yes
                                            
