インストールの仕方
================================================
[`Installation <https://www.rdkit.org/docs/Install.html#installation>`__]

以下は様々なインストール方法を多数示していますが、複雑さの度合いが異なります。

anaconda pythonを使ったクロスプラットフォームな方法(最速のインストール方法)
------------------------------------------------------------------------------------------------------------
[`Cross-platform under anaconda python(fastest install) <https://www.rdkit.org/docs/Install.html#cross-platform-under-anaconda-python-fastest-install>`__]

anacondaの紹介
------------------------------------------------------------------------------------------------------------
[`Introduction to anaconda <https://www.rdkit.org/docs/Install.html#introduction-to-anaconda>`__]

condaはオープンソースで、クロスプラットフォームのソフトウェアパッケージマネージャーです。
ソフトウェアコンポーネントのパッケージ化と配布をサポートしており、独立した実行環境内部でのインストールを管理します。
pipやvirtaulenvといくつかの類似性がありますが、より“Pythonに依存しない”形になるようデザインされていて、バイナリパッケージの配布と依存関係により適しています。

condaを手にいれる方法
------------------------------------------------------------------------------------------------------------
[`How to get conda <https://www.rdkit.org/docs/Install.html#how-to-get-conda>`__]

condaを手に入れる最も簡単な方法は\ `Anaconda Python distribution <https://conda.io/docs/user-guide/install/index.html>`__\ の一部分としてインストールすることです。
その他の方法として（少し使用方法が複雑になりますが）より小さく、必要なものに絞られた\ `Miniconda <https://conda.io/miniconda.html>`__\ をインストールするという選択肢もあります。
condaのソースコードレポジトリは\ `github <https://github.com/conda>`__\ で手に入ります。
追加のドキュメントもプロジェクトの\ `website <https://conda.io/docs/>`__\ 上で提供されています.

Condaを使ってRDKitをインストールする方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`How to install RDKit with Conda <https://www.rdkit.org/docs/Install.html#how-to-install-rdkit-with-conda>`__]

RDKitがインストールされた新しいconda環境を作るには、次のようなコマンドを一つ打ち込むだけでOKです:

.. code:: console

   $ conda create -c rdkit -n my-rdkit-env rdkit

最後に、同じシェル上で相当するPythonインタープリターが利用可能になるように新しい環境をアクティベートする必要があります:

.. code:: console

   $ source activate my-rdkit-env

もし、何らかの理由でうまくいかなかった場合は次を試してください:

.. code:: console

   $ cd [anaconda folder]/bin
   $ source activate my-rdkit-env

Windowsユーザーは少し違うコマンドが必要です:

.. code:: console

   C:\> activate my-rdkit-env

conda-forgeパッケージ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`conda-forge package <https://www.rdkit.org/docs/Install.html#conda-forge-package>`__]

\ `conda-forge <https://conda-forge.org/#about>`__\ のRDKitパッケージも利用可能です。
すでに他のconda-forgeパッケージを使っているユーザーにはもしかしたらより簡単なインストール方法かもしれません。
このパッケージは次のようにインストールすることができます:

.. code:: console

   $ conda install -c conda-forge rdkit

Condaを使ってソースからビルドする方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`How to build from source with Conda <https://www.rdkit.org/docs/Install.html#how-to-build-from-source-with-conda>`__]

condaを使ってソースからビルドする方法の詳細は、\ `conda-rdkit repository <https://github.com/rdkit/conda-rdkit>`__\ を見てください。

macOS 10.12 (Sierra): Python 3環境
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`macOS 10.12 (Sierra): Python 3 environment <https://www.rdkit.org/docs/Install.html#macos-10-12-sierra-python-3-environment>`__]

次のコマンドでmacOS Sierra と Python3のための開発環境を作ることができます。
\ `conda <http://conda.pydata.org/miniconda.html>`__\ からMiniconda3-latest-MacOSX-x86_64.shをダウンロードし次のコマンドを実行してください:

.. code:: console

   bash Miniconda3-latest-MacOSX-x86_64.sh
   conda install numpy matplotlib
   conda install cmake cairo pillow eigen pkg-config
   conda install boost-cpp boost py-boost

オプションとして、便利な開発用ツールとして、次のパッケージを付け足すこともできます:

.. code:: console

   pip install yapf==0.11.1
   pip install coverage==3.7.1

そして、通常のビルドの方法を実施します。PYTHON_INCLUDE_DIRをcmakeコマンドにセットする必要があります。

.. code:: console

   PYROOT=<path to miniconda3>
   cmake -DPYTHON_INCLUDE_DIR=$PYROOT/include/python3.6m  \
     -DRDK_BUILD_AVALON_SUPPORT=ON \
     -DRDK_BUILD_CAIRO_SUPPORT=ON \
     -DRDK_BUILD_INCHI_SUPPORT=ON \
     ..

\ ``make``\ と\ ``make install``\ が最後までうまくいったら、次のコマンドでテストを行ってください:

.. code:: console

   RDBASE=$RDBASE DYLD_FALLBACK_LIBRARY_PATH="$RDBASE/lib:$PYROOT/lib" PYTHONPATH=$RDBASE ctest

これは最近のmacOSバージョンに導入された\ `System Integrity Protection SIP <https://en.wikipedia.org/wiki/System_Integrity_Protection>`__\ のため、必要となります。

Linux x86_64: Python 3環境
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Linux x86_64: Python 3 environment <https://www.rdkit.org/docs/Install.html#linux-x86-64-python-3-environment>`__]

次のコマンドでLinux x86_64 と Python3のための開発環境を作ることができます。

まずは\ `Anaconda <https://www.anaconda.com/download/#linux>`__\ から最新のanacondaインストーラーをダウンロードしてください。
そして、必要なパッケージをインストールしてください:

.. code:: console

   bash Anaconda3-5.2.0-x86_64.sh
   conda install -y cmake cairo pillow eigen pkg-config
   conda install -y boost-cpp boost py-boost

Numpyとmatplotlibはanacondaの基本インストールの一部としてすでに含まれています。
最新のboostライブラリが、現在、anacondaのデフォルトよりも新しいGLIBICバージョンを使ってビルドされているため、
より新しいバージョンにアップデートする必要があります:

.. code:: console

   conda install -y gxx_linux-64

この時点で、ビルドしたい場所にRDKitレポジトリをクローンし、ビルドを開始できる必要があります。anacondaがNumpyパッケージの中に隠してしまうので、
RDKitがnumpyのヘッダを見つけられるようにパスを示す必要があります:

.. code:: console

   git clone https://github.com/rdkit/rdkit.git
   cd rdkit
   mkdir build && cd build
   cmake .. -DPy_ENABLE_SHARED=1 \
       -DRDK_INSTALL_INTREE=ON \
       -DRDK_INSTALL_STATIC_LIBS=OFF \
       -DRDK_BUILD_CPP_TESTS=ON \
       -DPYTHON_NUMPY_INCLUDE_PATH="$CONDA_PREFIX/lib/python3.6/site-packages/numpy/core/include"

最後に、\ ``make``\ 、\ ``make install``\ 、\ ``ctest``\ を実行してください。

conda環境からPostgreSQLとRDKit PostgreSQLカートリッジを使ってインストールする方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Installing and using PostgreSQL and the RDKit PostgreSQL cartridge from a conda environment <https://www.rdkit.org/docs/Install.html#installing-and-using-postgresql-and-the-rdkit-postgresql-cartridge-from-a-conda-environment>`__]

conda pythonディストリビューションがシステムのPythonのバージョンと異なるため、
PostgreSQLとPostgreSQL pythonクライアントはcondaを介してインストールするのが最も簡単です。

環境をアクティベートした状態で、次のコマンドを実行するだけです:

.. code:: console

   conda install -c rdkit rdkit-postgresql

condaパッケージのPostgreSQLバージョンは\ ``[conda folder]/envs/my-rdkit-env/bin``\ の中にある
initdbコマンドを走らせて初期化する必要があります：

.. code:: console

   [conda folder]/envs/my-rdkit-env/bin/initdb -D /folder/where/data/should/be/stored

初期化後、ターミナルから次のコマンドでPostgreSQLを実行することができます:

.. code:: console

   [conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored

多くの場合、代わりにdaemonとしてPostgreSQLを走らせる必要があると思います。これを行う一つの方法はsupervisorを使うことです。
より詳細な情報とsupervisorのインストール方法は\ `ここ <http://supervisord.org/>`__\ で手に入ります。
必要なコンフィギュレーションファイルは以下ようなものを見てください:

.. code:: console

   [program:postgresql]
   command=[conda folder]/envs/my-rdkit-env/bin/postgres -D /folder/where/data/should/be/stored
   user=[your username]
   autorestart=true

PostgreSQLが起動し、走り始めたら、conda環境をアクティベートすれば、通常のPostgreSQLコマンドは全て実行可能となります。
したがってデータベースを作るには次を実行してください:

.. code:: console

   createdb my_rdkit_db
   psql my_rdkit_db
   # create extension rdkit;

もし異なる環境でPostgreSQLを複数インストールして使おうとしているなら、\ `PostgreSQL configファイルを編集 <https://opensourcedbms.com/dbms/running-multiple-postgresql-9-2-instances-on-one-server-in-centos-6rhel-6fedora/>`__\ して、
異なるpidファイルとunixソケット、およびポートをセットアップする必要があります。
上記のコンフィギュレーションの場合、これらのファイルは\ ``/folder/where/data/should/be/stored``\ で見つけることができます。

Linux and OS X
-------------------------------------------------
[`Linux and OS X <https://www.rdkit.org/docs/Install.html#linux-and-os-x>`__]

レポジトリからインストールする方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Installation from repositories <https://www.rdkit.org/docs/Install.html#installation-from-repositories>`__]

Ubuntu 12.04 とそれ以降のバージョン
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Ubuntu 12.04 and later <https://www.rdkit.org/docs/Install.html#ubuntu-12-04-and-later>`__]

Debichemチームの努力のおかげで、RDKitはUbuntuレポジトリを介して手に入れることができます。インストールするには次を実行してください:

.. code:: console

   sudo apt-get install python-rdkit librdkit1 rdkit-data

Fedora、CentOSとRHEL
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Fedora, CentOS, and RHEL <https://www.rdkit.org/docs/Install.html#fedora-centos-and-rhel>`__]

Gianluca Sfornaのおかげで、RDKitのためのバイナリRPMが公式のFedoraレポジトリの一部となっています:
https://admin.fedoraproject.org/pkgdb/package/rpms/rdkit/

OS X
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`OS X <https://www.rdkit.org/docs/Install.html#os-x>`__]

Eddie CaoがRDKitを簡単にビルドするのに使うことができるhomebrewのフォーミュラを作りました:
https://github.com/rdkit/homebrew-rdkit

ソースからビルドする方法
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Building from Source <https://www.rdkit.org/docs/Install.html#building-from-source>`__]

\ ``2018_03``\ リリースから、RDKitのコアとなるC++コードはモダンなC++で書かれています。このリリースの場合C++11です。
つまり、ビルドするのに使うコパイラとして古いものを使うことが完全にできなくなったということです。最低限のテスト済みのバージョンは次のものです:

-  g++ v4.8: SLNパーサーコードはv4.8でビルドすることができないことに注意してください。古いコンパイラが使われると自動的に無効にされます。
-  clang v3.9: もっと古いバージョンも機能するかもしれませんが試していません。
-  Visual Studio 2015: もっと古いバージョンも機能するかもしれませんが試していません。

ソースから前もって必要となるものをインストールする方法
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Installing prerequisites from source <https://www.rdkit.org/docs/Install.html#installing-prerequisites-from-source>`__]

-  必要なパッケージ:
-  cmake。version
   3.1（以上）が必要です。Linuxディストリビューションが適したパーッケージを持っていなければhttp://www.cmake.org
   を参照してください。
-  Pythonラッパーを使おうと思っているなら次が必要です

   -  Pythonヘッダ。おそらくLinuxディストリビューションにpython-devパッケージ(あるいはそれが必要とするものを何でも）をインストールする必要があることを意味します。
   -  sqlite3。共有ライブラリも必要です。sqlite3-devパッケージをインストールする必要が有るかもしれません。
   -  Numpy
      (http://www.scipy.org/NumPy)をインストールしている必要があります。

         | **注意**
         |
         | OS XでXCode4を使ってビルドするには、XCode4に付随するnumpyのバージョンに問題があるようです。応急処置の方法は下のセクション(FAQセクション)を参照してください。

Boostのインストール
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
[`Installing Boost <https://www.rdkit.org/docs/Install.html#installing-boost>`__]

linxusディストリビューションに、Pythonとシリアライズ化ライブラリを含む、バージョン1.58以上のboost-develパッケージがあれば次のステップは必要ありません。

  | **注意**
  |
  | プレインストールされたboostライブラリのバージョンがある状態で、自分のバージョンを使いたい場合はコードをビルドする際に気をつけてください。我々は少なくとも一つ、Fedoraシステム上で、ユーザーがインストールしたバージョンのboostを使ってcmakeでコンパイルし、システムのバージョンにリンクさせた、という事例を知っていますが、この時はセグメンテーションフォールトになりました。応急処置の方法は下のセクション(FAQセクション)を参照してください。

-  \ `boostのウェブサイト <http://www.boost.org/>`__\ からboostのソースディストリビューションをダウンロード
-  マシンからソースを取得する（例えば\ ``/usr/local/src/boost_1_58_0``\ ）
-  必要なboostライブラリをビルドする。boostのサイトに\ `詳細な説明 <http://www.boost.org/doc/libs/1_58_0/more/getting_started/index.html>`__\ がありますが、全体の流れは次のようになります：
-  \ ``cd \$BOOST``\
-  Pythonラッパーを使いたい場合は:\ ``./bootstrap.sh --with-libraries=python,serialization``\
-  Pythonラッパーを使わない場合は:\ ``./bootstrap.sh --with-libraries=serialization``\
-  \ ``./b2 install``\

このステップで問題が生じた場合はboostの\ `インストールマニュアル <http://www.boost.org/more/getting_started/unix-variants.html>`__\ を参照してください。

Building the RDKit
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Building the RDKit <https://www.rdkit.org/docs/Install.html#building-the-rdkit>`__]

ソースを取得します。ここではtar.gzを使っていますがgitを使うこともできます。

.. code:: console

   wget https://github.com/rdkit/rdkit/archive/Release_XXXX_XX_X.tar.gz

-  あらかじめ必要なものがインストールされているか確認しておいてください
-  環境変数:

   -  RDBASE: RDKitディストリビューションのルートディレクトリ(例 ~/RDKit)
   -  *Linux:* LD_LIBRARY_PATH:
      $RDBASE/libとboost共有ライブラリがインストールされている場所が確実に含まれるようにしてください
   -  *OS X:* DYLD_LIBRARY_PATH:
      $RDBASE/libとboost共有ライブラリがインストールされている場所が確実に含まれるようにしてください
   -  Pythonラッパーを使おうと思っているなら次が必要です:
      -  PYTHONPATH: $RDBASEを含むことを確認してください

-  ビルド:
-  $RDBASE に cd
-  \ ``mkdir build``\
-  \ ``cd build``\
-  \ ``cmake ..``\ : デフォルトのバージョンではないPythonを指定する必要がある場合か、標準的ではない場所にboostがある場合は、ビルドを設定する際に下のセクションを参照してください
-  \ ``make``\ :（デフォルトでは）これで全てのライブラリとリグレッションテスト、ラッパーをビルドされます。
-  \ ``make install``\

FAQと解決方法のリストは下を見てください。

ビルドの検証 (オプションですが、推奨します)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Testing the build (optional, but recommended) <https://www.rdkit.org/docs/Install.html#testing-the-build-optional-but-recommended>`__]

-  $RDBASE/build にcd、次に\ ``ctest``\ を実行
-  これでおしまい！

発展的内容
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Advanced <https://www.rdkit.org/docs/Install.html#advanced>`__]

Boostのインストールの代わりとなるものを明示する
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
[`Specifying an alternate Boost installation <https://www.rdkit.org/docs/Install.html#specifying-an-alternate-boost-installation>`__]

cmakeにboostライブラリとヘッダファイルの場所を教える必要があります:
boostを\ ``/opt/local``\ においている場合、cmakeの呼び出しは次のようになります:

.. code:: console

   cmake -DBOOST_ROOT=/opt/local ..

システムがインストールしたboostがあるシステム上で、自分のboostを使っている場合は注意してください。
通常、cmakeコマンドに引数\ ``-D Boost_NO_SYSTEM_PATHS=ON``\ も含めた方が良いです。

Pythonのインストールの代わりとなるものを明示する
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
[`Specifying an alternate Python installation <https://www.rdkit.org/docs/Install.html#specifying-an-alternate-python-installation>`__]

デフォルトでインストールされたPythonを使っていない場合は、リンクすべきPythonライブラリの場所とPythonヘッダファイルの場所をcmakeに教える必要があります。

これがサンプルコマンドラインです:

.. code:: console

   cmake -D PYTHON_LIBRARY=/usr/lib/python2.7/config/libpython2.7.a -D PYTHON_INCLUDE_DIR=/usr/include/python2.7/ -D PYTHON_EXECUTABLE=/usr/bin/python ..

正しいPythonが、あなたのPATHの最初のバージョンなら、\ ``PYTHON_EXECUTABLE``\ の部分はオプションです。

Pythonのラッパーを無効化する
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
[`Disabling the Python wrappers <https://www.rdkit.org/docs/Install.html#disabling-the-python-wrappers>`__]

Pythonラッパーのビルドを完全に無効化することができます:

.. code:: console

   cmake -DRDK_BUILD_PYTHON_WRAPPERS=OFF ..

オススメの追加項目
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
[`Recommended extras <https://www.rdkit.org/docs/Install.html#recommended-extras>`__]

-  cmakeコマンドラインに\ ``-DRDK_BUILD_INCHI_SUPPORT=ON``\ 引数を加えることでInChI文字列とInChIキーを生成するためのサポートを有効にすることができます。
-  cmakeコマンドラインに\ ``-DRDK_BUILD_AVALON_SUPPORT=ON``\ 引数を加えることで
   Avalonツールキットのサポートを有効にすることができます。
-  構造描画のために高画質のPNGを生成できるようにしたいなら、cairoをシステムにインストールし、cairoのサポートを有効にしてRDKitをビルドする必要があります:\ ``-DRDK_BUILD_CAIRO_SUPPORT=ON``\
-  3D記述子を使えるようにしたいなら、eigen3のコピーがインストールされている状態にする必要があります。ほとんどのOSは適切なパッケージを持っています。

JAVAラッパーのビルド
''''''''''''''''''''''''''''''''''''''''''''''''''
[`Building the Java wrappers <https://www.rdkit.org/docs/Install.html#building-the-java-wrappers>`__]

*追加で必要となるもの*

-  SWIG >v2.0: http://www.swig.org

*ビルドの方法*

-  cmakeを呼び出すときに\ ``-D RDK_BUILD_SWIG_WRAPPERS=ON``\ を引数に加えてください。例えば:\ ``cmake -D RDK_BUILD_SWIG_WRAPPERS=ON ..``\
-  ビルドとインストールは通常makeを使います。
   \ ``$RDBASE/Code/JavaWrappers/gmwrapper``\ ディレクトリに3つの必要なファイルが入っています:
   libGraphMolWrap.so (libGraphMolWrap.jnilib on OS X)とorg.RDKit.jar、そしてorg.RDKitDoc.jarです。

*ラッパーを使う方法*

ラッパーを使うには3つのファイルが同じディレクトリに入っている必要があり、
そしてそれがCLASSPATHとjava.library.pathに含まれている必要があります。Jythonを使った例を示します:

.. code:: console

   % CLASSPATH=$CLASSPATH:$RDBASE/Code/JavaWrappers/gmwrapper/org.RDKit.jar; jython -Djava.library.path=$RDBASE/Code/JavaWrappers/gmwrapper
   Jython 2.2.1 on java1.6.0_20
   Type "copyright", "credits" or "license" for more information.
   >>> from org.RDKit import *
   >>> from java import lang
   >>> lang.System.loadLibrary('GraphMolWrap')
   >>> m = RWMol.MolFromSmiles('c1ccccc1')
   >>> m.getNumAtoms()
   6L

オプションのパッケージ
'''''''''''''''''''''''''''''''''''''''''''''''''''''
[`optional packages <https://www.rdkit.org/docs/Install.html#optional-packages>`__]

-  RDKit
   InChIのサポートをインストールしたいなら\ ``$RDBASE/External/INCHI-API/README``\ のインストラクションに従ってください。
-  RDKit
   Avalonツールキットのサポートをインストールしたいなら\ ``$RDBASE/External/AvalonTool/README``\ のインストラクションに従ってください。
-  PostgreSQLカートリッジをビルドしインストールしたいなら\ ``$RDBASE/Code/PgSQL/rdkit/README``\ のインストラクションに従ってください。

よくある問題
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Frequently Encountered Problems <https://www.rdkit.org/docs/Install.html#frequently-encountered-problems>`__]

以下の各事例ではパスの特定の部分を\ ``...``\ で置き換えています。

*問題:*

.. code:: console

   Linking CXX shared library libSLNParse.so
   /usr/bin/ld: .../libboost_regex.a(cpp_regex_traits.o): relocation R_X86_64_32S against `std::basic_string<char, std::char_traits<char>, std::allocator<char> >::_Rep::_S_empty_rep_storage' can not be used when making a shared object; recompile with -fPIC
   .../libboost_regex.a: could not read symbols: Bad value
   collect2: ld returned 1 exit status
   make[2]: *** [Code/GraphMol/SLNParse/libSLNParse.so] Error 1
   make[1]: *** [Code/GraphMol/SLNParse/CMakeFiles/SLNParse.dir/all] Error 2
   make: *** [all] Error 2

*解決方法:*

cmakeを呼び出すときに次を引数に加えてください:\ ``-DBoost_USE_STATIC_LIBS=OFF``\

さらに情報が欲しい場合はこちら：http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01119.html

*問題:*

.. code:: console

   .../Code/GraphMol/Wrap/EditableMol.cpp:114:   instantiated from here
   .../boost/type_traits/detail/cv_traits_impl.hpp:37: internal compiler error: in make_rtl_for_nonlocal_decl, at cp/decl.c:5067

   Please submit a full bug report, with preprocessed source if appropriate. See \<URL:<http://bugzilla.redhat.com/bugzilla>\> for instructions. Preprocessed source stored into /tmp/ccgSaXge.out file, please attach this to your bugreport. make[2]: **\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/EditableMol.cpp.o] Error 1 make[1]:**\* [Code/GraphMol/Wrap/CMakeFiles/rdchem.dir/all] Error 2 make: *\** [all] Error 2

*解決方法:*

\ ``Code/GraphMol/Wrap/EditableMol.cpp``\ の一番上に\ ``#define BOOST_PYTHON_NO_PY_SIGNATURES``\ を加えてください。

さらに情報が欲しい場合はこちら:
http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01178.html

*問題:*

システムには\ ``/usr/lib``\ にインストールされたバージョンのboostがあるが、RDKitにはより最新のものを使わせたい。

*解決方法:*

\ ``-D Boost_NO_SYSTEM_PATHS=ON``\ 引数を渡すことで解決できます。

*問題:*

OS XでXCode4を使ってビルドする場合

XCode4とともに配布されているNumpyのバージョンによって問題が起きているように見えるので、新しいコピーをビルドする必要があります。

*解決方法:*

Numpyのコピーを手に入れて、次のようにroot:をroot:としてビルドしてください:

.. code:: console

   export MACOSX_DEPLOYMENT_TARGET=10.6
   export LDFLAGS="-Wall -undefined dynamic_lookup -bundle -arch x86_64"
   export CFLAGS="-arch x86_64"
   ln -s /usr/bin/gcc /usr/bin/gcc-4.2
   ln -s /usr/bin/g++ /usr/bin/g++-4.2
   python setup.py build
   python setup.py install

確実に新しいNumpyがビルドで使われるようにしてしてください:

.. code:: console

   PYTHON_NUMPY_INCLUDE_PATH /Library/Python/2.6/site-packages/numpy/core/include

また、PYTHONPATHの最初になるようにしてください:

.. code:: console

   export PYTHONPATH="/Library/Python/2.6/site-packages:$PYTHONPATH"

これで、安全にboostとRDKitをビルドすることができます。

Windows
---------------------------------------------------------------------
[`Windows <https://www.rdkit.org/docs/Install.html#windows>`__]

必要条件
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Prerequisites <https://www.rdkit.org/docs/Install.html#prerequisites>`__]

-  3.6+ (http://www.python.org/)
-  numpy （http://numpy.scipy.org/ あるいは\ ``pip install numpy``\ ）。
   win64用のバイナリはここで手に入る: http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy
-  Pillow:（https://python-pillow.github.io/> あるいは\ ``pip install Pillow``\ ）

オススメの追加項目
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Recommeded extras <https://www.rdkit.org/docs/Install.html#recommended-extras>`__]

-  aggdraw:Pythonで高画質の描画を行うためのライブラリ。ダウンロードの方法はこちら:
   http://effbot.org/zone/aggdraw-index.htm
   新しい描画コード（2008年5月時点)はaggdraw v1.2a3 でテストした。alpha labelにも関わらず、コードは安定していてかつ機能する。
-  matplotlib: Pythonから科学的なプロットを作成するためのライブラリ。
   http://matplotlib.sourceforge.net/
-  ipython :
   Python（とそれ以上に多くのこと）のためのとても役に立つインタラクティブなシェル
   http://ipython.scipy.org/dist/

   -  win32all: PythonのためのWindowsの拡張
      http://sourceforge.net/projects/pywin32/

RDKitバイナリのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Installation of RDKit binaries <https://www.rdkit.org/docs/Install.html#installation-of-rdkit-binaries>`__]

-  適切なWindows用のバイナリビルドを
   https://github.com/rdkit/rdkit/releases から入手する
-  zipファイルをどこかから取り出す(i.e. \ ``C:\``\ )。名前にスペースを入れないこと。
-  以下は\ ``C:\RDKit_2015_09_2``\ にインストールされていると仮定して進める
-  次の環境変数を設定する：
-  RDBASE:\ ``C:\RDKit_2015_09_2``\
-  PYTHONPATH:\ ``%RDBASE%``\
   もしすでにPYTHONPATHがあるなら、末尾に\ ``;%RDBASE%``\ を付け足す。。
-  PATH:\ ``;%RDBASE%\lib``\ を末尾に付け足す。

Win7システムではDLLが見つからないことでトラブルに見舞われるかもしれません。
メーリングリストのスレッド(http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg01632.html)
を参照してください。必要なDLLはここからダウンロードできます: http://www.microsoft.com/en-us/download/details.aspx?id=5555

ソースからのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[`Installation from source <https://www.rdkit.org/docs/Install.html#installation-from-source>`__]

追加のソフトウェアのインストール
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Extra software to install <https://www.rdkit.org/docs/Install.html#extra-software-to-install>`__]

-  Microsoft Visual C++ :
   コミュニティバージョンに必要なものは全て揃っていて、無料でダウンロードすることができます(https://www.visualstudio.com/vs/community)。
   大きなインストールになるので時間がかかります。RDKitはVisual Studio 2015と2017を使ってビルドされています、より最新のバージョンの方が良いでしょう。
-  cmake : (http://www.cmake.org/cmake/resources/software.html)
   インストールする必要があります。
-  boost : boostライブラリのコンパイル済みのバージョンを http://sourceforge.net/projects/boost/files/boost-binaries/
   からダウンロートして使うことを強くお勧めします。インストーラーを実行するときに必要なバイナリライブラリはPythonとserializationだけです。
   ソースからboostをインストールしたいなら、http://www.boost.org からコピーをダウンロードして、
   ドキュメントの“Getting Started”セクションの指示に従ってください。
   確実にライブラリとヘッダが\ ``C:\boost``\ にインストールされるようにしてください。
-  a git client : *これはRDKitのdevelopment versionでビルドしようと思っている時だけ必要になります。*
   http://git-scm.com/downloads からダウンロードすることができます。
   gitはMicrosoft Visual Studio2015のオプションのアドオンとしても含まれています。

セットアップと準備
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Setup and Preparation <https://www.rdkit.org/docs/Install.html#setup-and-preparation>`__]

このセクションはPythonが\ ``C:\Python27``\ に、boostライブラリが\ ``C:\boost``\ にインストールされていて、
RDKitを\ ``C:\RDKit``\ という名前のディレクトリからビルドしようとしていると仮定して進めます。
以上の条件のどれか一つでも違っているなら、対応するパスを変更してください。

-  名前にスペースを含むパスにインストールしようとしているなら、確実に環境変数の定義に適切に引用を使うようにしてください。
-  RDKitのdevelopment
   versionを使おうとしているなら、gitを使っている現在のRDKitのソースのコピーを手に入れてください。
   コマンドラインクライアントを使っているなら、コマンドは\ ``git clone https://github.com/rdkit/rdkit.git C:\RDKit``\ です。
-  リリースバージョンのRDKitを使おうとしているなら、最新のリリースのコピーを手に入れ、\ ``C:\RDKit``\ ディレクトリに取り出してください。
-  必要な環境変数を設定します:
-  \ ``RDBASE = C:\RDKit``\
-  \ ``C:\Python27``\ がPATHに入っていることを確認してください
-  \ ``C:\RDKit\lib``\ がPATHに入っていることを確認してください
-  \ ``C:\boost\lib``\  がPATHに入っていることを確認してください
-  \ ``C:\RDKit``\ がPYTHONPATHに入っていることを確認してください

コマンドラインからのビルド（推奨の方法）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Building from the command line (recommended) <https://www.rdkit.org/docs/Install.html#building-from-the-command-line-recommended>`__]

-  \ ``C:\RDKit\build``\ ディレクトリを作ってそこにcdしてください
-  cmakeを実行してください。InChIとAvalonツールキットのソースをそれぞれInChI
   TrustとSourceForgeレポジトリからダウンロードし、インストール済みのバージョンのPostgreSQLのためのPostgreSQL
   カートリッジをビルドする、64bit版Windowsのための基本的なコマンドラインの例は次のようになります。
   \ ``cmake -DRDK_BUILD_PYTHON_WRAPPERS=ON -DBOOST_ROOT=C:/boost -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_AVALON_SUPPORT=ON -DRDK_BUILD_PGSQL=ON -DPostgreSQL_ROOT="C:\Program Files\PostgreSQL\9.5" -G"Visual Studio 14 2015 Win64" ..``\
-  コードをビルドしてください。コマンドラインの例は次のようになります:\ ``C:/Windows/Microsoft.NET/Framework64/v4.0.30319/MSBuild.exe /m:4 /p:Configuration=Release INSTALL.vcxproj``\
-  PostgreSQLサポートでビルドした場合、管理者権限でシェルを開き、PostgreSQLサービスをとめ、インストール用のスクリプト\ ``pgsql_install.bat``\ を実行し、そしてPostgreSQLサービスを再起動する必要があります。
   （より詳細には\ ``%RDBASE%\Code\PgSQL\rdkit\README``\ を参照してください）:

   -  \ ``"C:\Program Files\PostgreSQL\9.5\bin\pg_ctl.exe" -N “postgresql-9.5” -D “C:\Program Files\PostgreSQL\9.5\data” -w stop``\
   -  \ ``C:\RDKit\build\Code\PgSQL\rdkit\pgsql_install.bat``\
   -  \ ``"C:\Program Files\PostgreSQL\9.5\bin\pg_ctl.exe" -N "postgresql-9.5" -D "C:\Program Files\PostgreSQL\9.5\data" -w start``\
   -  PostgreSQLサービスを再起動する前に、RDKitをビルドする際に使ったBoostライブラリがシステムのPATHにあることを確認してください。さもないとPostgreSQLは\ ``rdkit``\ の拡張を作るのに失敗し、次のような紛らわしいエラーメッセージを出します:\ ``ERROR: could not load library "C:/Program Files/PostgreSQL/9.5/lib/rdkit.dll": The specified module could not be found.``\

ビルドの検証（オプションですが、推奨します）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
[`Testing the Build (optional, but recommended) <https://www.rdkit.org/docs/Install.html#testing-the-build-optional-but-recommended>`__]

-  \ ``C:\RDKit\build``\ にcdし、ctestを実行してください。PostgreSQLサポート内部でビルドした場合は、現在のログインユーザーが、データベースの作成とスーパーユーザーの権限を持ったPostgreSQLユーザーである必要があることに注意してください。さもないとPostgreSQLのテストは失敗します。認証のための便利なオプションは、ctestを実行しているシェルの現在のログインユーザーのPostgreSQLパスワードに、\ ``PGPASSWORD``\ 環境変数をセットすることです。
-  これでおしまい

ライセンス
--------------------------------------------------------
[`License <https://www.rdkit.org/docs/Install.html#license>`__]

この文書の著作権は copyright (C) 2013-2018 by Greg Landrum に所属しています。

この文書はCreative Commons Attribution-ShareAlike 4.0 Licenseのもとでライセンスされています。
このライセンスを見るためには http://creativecommons.org/licenses/by-sa/4.0/ にアクセスするか、
Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.　に手紙を送ってください。

このライセンスの意図はRDKitそのものの意図と似ています。
簡単に言えば“これを使ってなんでもやりたいことをやっていいですが、私たちの功績についても言及してください”
