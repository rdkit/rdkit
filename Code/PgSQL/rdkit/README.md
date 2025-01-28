PgSQL/rdkit is a PostgreSQL contribution module, which implements
*mol* datatype to describe molecules and *fp* datatype for fingerprints,
basic comparison operations, similarity operation (tanimoto, dice) and  index
support (using GiST indexing framework).

**Compatibility**: PostgreSQL 9.1+
If using PostgreSQL packages from your linux distribution, be sure
that you have the `-devel` package installed.

**Installation**:
To build the optional PostgreSQL cartridge, set the CMake switch
`RDK_BUILD_PGSQL` to `ON`:
```
    -D RDK_BUILD_PGSQL=ON
```

Set `RDK_PGSQL_STATIC` to `ON` to statically link all RDKit libraries; to `OFF`
to dynamically link them:
```
    -D RDK_PGSQL_STATIC=OFF
```

# Linux

## Cartridge installation into conda PostgreSQL
If you build the PostgreSQL cartridge against the conda `postgresql` package,
CMake will find it without having to specify its location through CMake
switches. Also installation will be straightforward as it will not require
root privileges.
Build RDKit in the usual way running `cmake`, followed by `make && make install`.
You should see the following message:
```
=====================================================================
PostgreSQL cartridge installation succeeded.
=====================================================================
Restart the PostgreSQL service before using the PostgreSQL
RDKit cartridge.
=====================================================================
```
To test the cartridge, first I initialize a database directory:
```
$ initdb /home/build/postgresdb
```
I edit the /home/build/postgresdb/postgresql.conf configuration file
to set the connection port to 6432 since on my system the default
`5432` is already in use, then I start the PostgreSQL server:
```
$ pg_ctl -D /home/build/postgresdb -l /home/build/postgresdb/logfile start
waiting for server to start.... done
server started
```
To test the RDKit cartridge, I create a a small test database:
```
$ createdb -p 6432 nci5K
```
And attempt to load the cartridge in the database:
```
$ psql -c 'create extension rdkit' -p 6432 -h /tmp nci5K
CREATE EXTENSION
```

## Cartridge installation into non-conda PostgreSQL
If PostgreSQL is installed in a location where CMake is unable to find it,
point the PostgreSQL-related CMake switches to the appropriate directories.
For example, these are the settings for PostgreSQL 15 on Ubuntu 20.04 LTS:
```
    -D RDK_BUILD_PGSQL=ON \
    -D PostgreSQL_CONFIG_DIR=/usr/lib/postgresql/15/bin \
    -D PostgreSQL_INCLUDE_DIR="/usr/include/postgresql" \
    -D PostgreSQL_TYPE_INCLUDE_DIR="/usr/include/postgresql/15/server" \
    -D PostgreSQL_LIBRARY="/usr/lib/x86_64-linux-gnu/libpq.so" \
```
Build RDKit in the usual way running `cmake`, followed by `make && make install`.
If PostgreSQL is installed in a directory that requires root privileges for
write access, the installation will fail:
```
=====================================================================
PostgreSQL cartridge installation failed:
=====================================================================
```
You will see a message explaining what you should do to successfully
install the cartridge. For example, on my WSL Ubuntu 20.04 LTS install I get:
```
=====================================================================
This might be due to insufficient privileges.
Check /home/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
for correctness of installation paths. If everything is OK, gain
administrator privileges, stop the PostgreSQL service, run
/home/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
to install the PostgreSQL RDKit cartridge, then start again
the PostgreSQL service.
=====================================================================
```
If you can gain sudo privileges, you can install the cartridge
as shown below.
```
# On SysV-based Linux distributions such as WSL2 Ubuntu 20.04
# replace with sudo service postgresql stop
$ sudo systemctl stop postgresql
 * Stopping PostgreSQL 15 database server
$ sudo sh /home/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
+ test -n
+ cp /home/build/build/rdkit/Code/PgSQL/rdkit/rdkit--4.2.sql /usr/share/postgresql/15/extension
+ cp /home/build/src/rdkit/Code/PgSQL/rdkit/rdkit.control /usr/share/postgresql/15/extension
+ cp /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--3.8--4.0.sql /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.0--4.0.1.sql /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.0--4.1.sql /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.0.1--4.1.sql /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.1--4.1.1.sql /home/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.1.1--4.2.sql /usr/share/postgresql/15/extension
+ cp /home/build/build/rdkit/Code/PgSQL/rdkit/librdkit.so /usr/lib/postgresql/15/lib/rdkit.so
# On SysV-based Linux distributions such as WSL2 Ubuntu 20.04
# replace with sudo service postgresql start
$ sudo systemctl start postgresql
 * Starting PostgreSQL 15 database server
```
To test the RDKit cartridge, I create a `build` user with superuser
privileges:
```
$ sudo -u postgres createuser -d -s -P build
Enter password for new role:
Enter it again:
```
Then I create a small test database:
```
$ createdb nci5K
```
And attempt to load the cartridge in the database:
```
$ psql -c 'create extension rdkit' nci5K
ERROR:  could not load library "/usr/lib/postgresql/15/lib/rdkit.so": /lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.29' not found (required by /home/build/install/rdkit_wsl64/lib/libRDKitAvalonLib.so.1)
```
The reason why I get the above error is that I built RDKit against
conda libraries, and the conda `libstdc++.so.6` is newer than Ubuntu's.
Therefore I need to make sure that PostgreSQL loads the cartridge
using the conda libraries preferentially, when possible.
For this, I need to adjust the PostgreSQL environment in
`/etc/postgresql/<version>/<cluster>/environment`.

Therefore, I set `LD_LIBRARY_PATH` in `/etc/postgresql/15/main/environment`
to include RDKit and conda libraries:
```
$ echo "LD_LIBRARY_PATH = '/home/build/install/rdkit_wsl64/lib:$CONDA_PREFIX/lib'" | sudo tee -a /etc/postgresql/15/main/environment
```
Then I restart PostgreSQL:
```
# On SysV-based Linux distributions such as WSL2 Ubuntu 20.04
# replace with sudo service postgresql restart
$ sudo systemctl restart postgresql
 * Restarting PostgreSQL 15 database server
```
As expected, this time the cartridge loads successfully:
```
$ psql -c 'create extension rdkit' nci5K
CREATE EXTENSION
```
The RDKit cartridge is now to ready to be used on the `nci5K` database.


# Windows

## Cartridge installation into conda PostgreSQL
If you build the PostgreSQL cartridge against the conda `postgresql` package,
CMake will find it without having to specify its location through CMake
switches. Also installation will be straightforward as it will not require
root privileges.
Build RDKit in the usual way: configure with `cmake`, then build and install:
```
>cmake --build . --parallel %NUMBER_OF_PROCESSORS% --target install --config Release
```
You should see the following message:
```
=====================================================================
PostgreSQL cartridge installation succeeded.
=====================================================================
Restart the PostgreSQL service before using the PostgreSQL
RDKit cartridge.
=====================================================================
```
To test the cartridge, first I initialize a database directory:
```
>initdb M:\postgresdb
```
I edit the M:\postgresdb\postgresql.conf configuration file
to set the connection port to 6433 since on my system the default
`5432` is already in use, then I start the PostgreSQL server:
```
>pg_ctl -D M:\postgresdb -l M:\postgresdb\logfile start
waiting for server to start.... done
server started
```
To test the RDKit cartridge, I create a a small test database:
```
>createdb -p 6433 nci5K
```
And attempt to load the cartridge in the database:
```
>psql -c 'create extension rdkit' -p 6433 nci5k
CREATE EXTENSION
```
Should you see an error message such as
```
ERROR:  could not load library "m:/a3/envs/rdkit_devel/Library/lib/rdkit.dll": The specified module could not be found.
```
this would indicate that some of the dependencies cannot be located in the
current `PATH`; for example, in my case I had to add the directory where
I chose to install RDKit DLLs as well as the MSYS2 directory containing
the `cairo-2.dll` dependency:
```
>set PATH=M:\install\rdkit_vs2022_conda_postgres\Lib;C:\msys64\mingw64\bin;%PATH%
```
Once you set your `PATH` appropriately, the cartridge should then load
without errors.

## Cartridge installation into non-conda PostgreSQL
if PostgreSQL is installed in a location where CMake is unable to find it,
point the PostgreSQL-related CMake switches to the appropriate directories.
For example, these are the settings for PostgreSQL 15 on Windows 11:
```
    -D RDK_INSTALL_DLLS_MSVC=ON ^
    -D RDK_BUILD_PGSQL=ON ^
    -D RDK_PGSQL_STATIC=OFF ^
    -D PostgreSQL_INCLUDE_DIR="C:/Program Files/PostgreSQL/15/include" ^
    -D PostgreSQL_TYPE_INCLUDE_DIR="C:/Program Files/PostgreSQL/15/include/server" ^
    -D PostgreSQL_LIBRARY="C:/Program Files/PostgreSQL/15/lib/postgres.lib" ^
```
Build RDKit in the usual way: configure with `cmake`, then build and install:
```
>cmake --build . --parallel %NUMBER_OF_PROCESSORS% --target install --config Release
```
If PostgreSQL is installed in a directory that requires root privileges for
write access, the installation will fail:
```
=====================================================================
PostgreSQL cartridge installation failed:
=====================================================================
```
You will see a message explaining what you should do to successfully
install the cartridge. For example, on my Windows 11 install I get:
```
=====================================================================
This might be due to insufficient privileges.
Check M:/build/rdkit_vs2022/Code/PgSQL/rdkit/pgsql_install.bat
for correctness of installation paths. If everything is OK, gain
administrator privileges, stop the PostgreSQL service, run
M:/build/rdkit_vs2022/Code/PgSQL/rdkit/pgsql_install.bat
to install the PostgreSQL RDKit cartridge, then start again
the PostgreSQL service.
=====================================================================
```
If you can open a CMD Window with administrator privileges, you can install
the cartridge as shown below.

Stop the PostgreSQL service:
```
>net stop postgresql-x64-15
The postgresql-x64-15 service is stopping.
The postgresql-x64-15 service was stopped successfully.
```
Install RDKit cartridge into PostgreSQL:
```
>M:/build/rdkit_vs2022/Code/PgSQL/rdkit/pgsql_install.bat

c:\Program Files>copy /Y "M:\build\rdkit_vs2022\Code\PgSQL\rdkit\rdkit--4.2.sql" "C:\PROGRA~1\POSTGR~1\15\share\extension"
        1 file(s) copied.

c:\Program Files>copy /Y "M:\src\rdkit\Code\PgSQL\rdkit\rdkit.control" "C:\PROGRA~1\POSTGR~1\15\share\extension"
        1 file(s) copied.

c:\Program Files>copy /Y "M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\*.sql" "C:\PROGRA~1\POSTGR~1\15\share\extension"
M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\rdkit--3.8--4.0.sql
M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\rdkit--4.0--4.0.1.sql
M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\rdkit--4.0.1--4.1.sql
M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\rdkit--4.1--4.1.1.sql
M:\build\rdkit_vs2022\Code\PgSQL\rdkit\update_sql\rdkit--4.1.1--4.2.sql
        5 file(s) copied.

c:\Program Files>copy /Y "M:\build\rdkit_vs2022\Code\PgSQL\rdkit\Release\rdkit.dll" "C:\PROGRA~1\POSTGR~1\15\lib\rdkit.dll"
        1 file(s) copied.
```
Make sure that the RDKit DLLs and their dependencies (in my case, conda
libraries and `cairo` libraries from MSYS2) are in the PostgreSQL `PATH`;
you can set this in the registry:
```
reg add HKLM\SYSTEM\CurrentControlSet\Services\postgresql-x64-15 /v Environment /t REG_MULTI_SZ /d "PATH=M:\install\rdkit_vs2022\Lib;M:\a3\envs\rdkit_devel\Library\bin;C:\msys64\mingw64\bin;%PATH%" /f
```
Start the PostgreSQL service:
```
>net start postgresql-x64-15
The postgresql-x64-15 service is starting.
The postgresql-x64-15 service was started successfully.
```
To test the RDKit cartridge, I create a `build` user with superuser
privileges:
```
>"C:\Program Files\PostgreSQL\15\bin\createuser.exe" -d -s -p 6432 -U postgres -P build
Enter password for new role:
Enter it again:
Password:
```
Then, in a normal CMD prompt as the build user, I create a small test database:
```
>"C:\Program Files\PostgreSQL\15\bin\createdb.exe" -p 6432 -W nci5K
Password:
```
And attempt to load the cartridge in the database:
```
"C:\Program Files\PostgreSQL\15\bin\psql.exe" -c "create extension rdkit" -p 6432 -W nci5K
Password:
CREATE EXTENSION
```
If you get an error such as
```
ERROR:  could not load library "C:/Program Files/PostgreSQL/15/lib/rdkit.dll": The specified module could not be found.
```
it means that one or more DLL that `rdkit.dll` depends on could not be found
in the path you set with the `reg add` command.
The best tools to troubleshoot missing DLL problems are
`DependenciesGui.exe` and `Procmon64.exe`; the latter in particular will show
exactly which DLL is trying to be loaded without success.


# macOS

## Cartridge installation into conda PostgreSQL
If you build the PostgreSQL cartridge against the conda `postgresql` package,
CMake will find it without having to specify its location through CMake
switches. Also installation will be straightforward as it will not require
root privileges.
Build RDKit in the usual way running `cmake`, followed by `make && make install`.
You should see the following message:
```
=====================================================================
PostgreSQL cartridge installation succeeded.
=====================================================================
Restart the PostgreSQL service before using the PostgreSQL
RDKit cartridge.
=====================================================================
```
To test the cartridge, first I initialize a database directory:
```
% initdb /Users/build/postgresdb
```
I edit the /Users/build/postgresdb/postgresql.conf configuration file
to set the connection port to 6432 since on my system the default
`5432` is already in use, then I start the PostgreSQL server:
```
% pg_ctl -D /Users/build/postgresdb -l /Users/build/postgresdb/logfile start
waiting for server to start.... done
server started
```
To test the RDKit cartridge, I create a a small test database:
```
% createdb -p 6432 nci5K
```
And attempt to load the cartridge in the database:
```
% psql -c 'create extension rdkit' -p 6432 nci5K
CREATE EXTENSION
```

## Cartridge installation into non-conda PostgreSQL
if PostgreSQL is installed in a location where CMake is unable to find it,
point the PostgreSQL-related CMake switches to the appropriate directories.
For example, these are the settings for PostgreSQL 15 on macOS 13.0.1:
```
    -D RDK_BUILD_PGSQL=ON \
    -D RDK_PGSQL_STATIC=OFF \
    -D PostgreSQL_INCLUDE_DIR=/Library/PostgreSQL/15/include \
    -D PostgreSQL_TYPE_INCLUDE_DIR=/Library/PostgreSQL/15/include/postgresql/server \
    -D PostgreSQL_LIBRARY=/Library/PostgreSQL/15/lib/libpq.dylib \
```
Build RDKit in the usual way running `cmake`, followed by `make && make install`.
If PostgreSQL is installed in a directory that requires root privileges for
write access, the installation will fail:
```
=====================================================================
PostgreSQL cartridge installation failed:
=====================================================================
```
You will see a message explaining what you should do to successfully
install the cartridge. For example, on my WSL Ubuntu 20.04 LTS install I get:
```
=====================================================================
This might be due to insufficient privileges.
Check /Users/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
for correctness of installation paths. If everything is OK, gain
administrator privileges, stop the PostgreSQL service, run
/Users/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
to install the PostgreSQL RDKit cartridge, then start again
the PostgreSQL service.
=====================================================================
```
If you can gain sudo privileges, you can install the cartridge
as shown below.
```
 sudo -u postgres /Library/PostgreSQL/15/bin/pg_ctl -D /Library/PostgreSQL/15/data stop
could not identify current directory: Permission denied
waiting for server to shut down.... done
server stopped
% sudo sh /Users/build/build/rdkit/Code/PgSQL/rdkit/pgsql_install.sh
+ test -n ''
+ cp /Users/build/build/rdkit/Code/PgSQL/rdkit/rdkit--4.2.sql /Library/PostgreSQL/15/share/postgresql/extension
+ cp /Users/build/src/rdkit/Code/PgSQL/rdkit/rdkit.control /Library/PostgreSQL/15/share/postgresql/extension
+ cp /Users/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--3.8--4.0.sql /Users/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.0--4.0.1.sql /Users/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.0.1--4.1.sql /Users/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.1--4.1.1.sql /Users/build/build/rdkit/Code/PgSQL/rdkit/update_sql/rdkit--4.1.1--4.2.sql /Library/PostgreSQL/15/share/postgresql/extension
+ cp /Users/build/build/rdkit/Code/PgSQL/rdkit/rdkit.so /Library/PostgreSQL/15/lib/postgresql/rdkit.so
% sudo -u postgres /Library/PostgreSQL/15/bin/pg_ctl -D /Library/PostgreSQL/15/data start
waiting for server to start....2023-03-19 21:45:13.230 CET [60023] LOG:  redirecting log output to logging collector process
2023-03-19 21:45:13.230 CET [60023] HINT:  Future log output will appear in directory "log".
 done
server started
```
To test the RDKit cartridge, I create a `build` user with superuser
privileges:
```
% sudo -u postgres /Library/PostgreSQL/15/bin/createuser -d -s -P build
Enter password for new role:
Enter it again:
Password:
```
Then I create a small test database:
```
% /Library/PostgreSQL/15/bin/createdb nci5K
Password:
```
And attempt to load the cartridge in the database:
```
% /Library/PostgreSQL/15/bin/psql -c 'create extension rdkit' nci5K
Password for user build:
ERROR:  could not load library "/Library/PostgreSQL/15/lib/postgresql/rdkit.so": dlopen(/Library/PostgreSQL/15/lib/postgresql/rdkit.so, 0x000A): Library not loaded: @rpath/libRDKitAvalonLib.1.dylib
  Referenced from: [...]
```
The reason why I get the above error is that some of the libraries
are in a location which is not listed among the `LC_RPATH`s of
`/Library/PostgreSQL/15/lib/postgresql/rdkit.so`, hence cannot be found:
```
% otool -l /Library/PostgreSQL/15/lib/postgresql/rdkit.so | grep LC_RPATH -A2
          cmd LC_RPATH
      cmdsize 56
         path /Users/build/a3/envs/rdkit_devel/lib (offset 12)
```
I can add another `LC_RPATH` using `install_name_tool`:
```
% sudo install_name_tool -add_rpath /Users/build/install/rdkit/lib /Library/PostgreSQL/15/lib/postgresql/rdkit.so
```
As expected, this time the cartridge loads successfully:
```
% /Library/PostgreSQL/15/bin/psql -c 'create extension rdkit' nci5K
Password for user build:
CREATE EXTENSION
```
The RDKit cartridge is now to ready to be used on the `nci5K` database.


**Further Reading**:

[Cartridge tutorial](https://rdkit.org/docs/Cartridge.html)
