%define version 1.95.5
%define release 1

Summary: Expat is an XML 1.0 parser written in C.
Name: expat
Version: %{version}
Release: %{release}
Copyright: MIT/X
Group: Utilities/parsers
URL: http://www.libexpat.org/
Source: http://download.sourceforge.net/expat/expat-%{version}.tar.gz
BuildRoot: /var/tmp/%{name}-buildroot

%description
Expat is an XML 1.0 parser written in C by James Clark.  It aims to be
fully conforming. It is not a validating XML parser.

%prep
%setup

%build
./configure
make lib xmlwf

%install
rm -rf $RPM_BUILD_ROOT
mkdir -p $RPM_BUILD_ROOT/usr/bin
mkdir -p $RPM_BUILD_ROOT/usr/lib
mkdir -p $RPM_BUILD_ROOT/usr/include
make install prefix=$RPM_BUILD_ROOT/usr
install -D xmlwf/xmlwf $RPM_BUILD_ROOT/usr/bin/xmlwf

%files
%doc COPYING Changes MANIFEST README doc/reference.html doc/style.css doc/valid-xhtml10.png
/usr/bin/xmlwf
/usr/lib
/usr/include/expat.h
/usr/man/man1/xmlwf.1.gz

%changelog
* Wed Sep 4 2002 Fred L. Drake, Jr. <fdrake@acm.org>
[Release 1.95.5-1]
- Updated for the 1.95.5 release.
- Updated URL for Expat home page to point to www.libexpat.org.
- Added "Valid XHTML 1.0" icon to the installed documentation.

* Sat Jun 29 2002 Fred L. Drake, Jr. <fdrake@acm.org>
[Release 1.95.4-1]
- Updated for the 1.95.4 release.

* Fri May 17 2002 Fred L. Drake, Jr. <fdrake@acm.org>
[Release 1.95.3-1]
- Updated for the 1.95.3 release.
- Added xmlwf man page to the list of files.

* Wed Jul 25 2001 Fred L. Drake, Jr. <fdrake@acm.org>
[Release 1.95.2-1]
- Updated for the 1.95.2 release.

* Sun Feb 18 2001 Sean Reifschneider <jafo-rpms@tummy.com>
[Release 1.95.1-1tummy]
- Updated to 1.95.1 release.
- Removed the "/usr/include/expat" directory for headers, as it now uses
  "expat.h" instead of "xmlparser.h".

* Thu Jan 25 2001 Sean Reifschneider <jafo-rpms@tummy.com>
[Release 1.1-3tummy]
- Moved xmlparse.h into "/usr/include/expat" directory to prevent conflict
  with w3c-libwww-devel package.

* Wed Sep 6 2000 Sean Reifschneider <jafo-rpms@tummy.com>
- Modified to install into /usr.
- Modified to use RPM_BUILD_ROOT instead of writing directly to install
  location.
