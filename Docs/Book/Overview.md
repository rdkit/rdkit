# An overview of the RDKit

## What is it?

### Open source toolkit for cheminformatics
-   Business-friendly BSD license
-   Core data structures and algorithms in C++
-   Python 3.x wrappers generated using Boost.Python
-   Java and C\# wrappers generated with SWIG
-   JavaScript wrappers of most-important functionality
-   2D and 3D molecular operations
-   Descriptor generation for machine learning
-   Molecular database cartridge for PostgreSQL
-   Cheminformatics nodes for KNIME (distributed from the KNIME community site: https://www.knime.com/rdkit)

### Operational:
- http://www.rdkit.org
- Supports Mac/Windows/Linux
- Major releases every 6 months, minor release about once a month
- Web presence:
    - Homepage: https://www.rdkit.org
      Documentation, links
    - Github (https://github.com/rdkit)
      Downloads, discussion, bug tracker, git repository
    - Sourceforge (https://sourceforge.net/projects/rdkit)
      Mailing lists
    - Blog (https://greglandrum.github.io/rdkit-blog/)
      Tips, tricks, random stuff
    - KNIME integration (https://github.com/rdkit/knime-rdkit)
      RDKit nodes for KNIME
- Mailing lists at https://sourceforge.net/p/rdkit/mailman/, searchable archives available for [rdkit-discuss](http://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/) and [rdkit-devel](http://www.mail-archive.com/rdkit-devel@lists.sourceforge.net/)
- Social media:
    - BlueSky: @rdkit.bsky.social
    - Mastodon: rdkit@mastodon.social
    - LinkedIn: https://www.linkedin.com/groups/8192558
    - Slack: https://rdkit.slack.com (invite required, contact Greg)

## Citing the RDKit

There is no official RDKit publication, our recommended citation is:
```
RDKit: Open-source cheminformatics. https://www.rdkit.org
```
We also recommend that you include the DOI for the version of the RDKit you used in the work. You can look these up here:
[https://doi.org/10.5281/zenodo.591637](https://doi.org/10.5281/zenodo.591637)


### Powered by RDKit
[![RDKit badge](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

If you use RDKit in one of your projects, you can show your support and help us track it by adding our badge.
Simply copy the code from one of the markup languages below and paste it in your README file:

<details>
  <summary>Markdown</summary>
  <div style="display: flex">
    <textarea rows="1" wrap="off" readonly style="display: inline-block; width: 100%; resize: none; overflow-y: hidden">[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    </textarea>
    <button onclick="copyBadgeMarkdown(this)" style="display: inline-block; padding: 0 8px">Copy</button>
  </div>
</details>

<details>
  <summary>reStructuredText</summary>
  <div style="display: flex">
    <textarea rows="1" wrap="off" readonly style="display: inline-block; width: 100%; resize: none; overflow-y: hidden">.. image:: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC
      :alt: Powered by RDKit
      :target: https://www.rdkit.org/
    </textarea>
    <button onclick="copyBadgeMarkdown(this)" style="display: inline-block; padding: 0 8px">Copy</button>
  </div>
</details>

<details>
  <summary>HTML</summary>
  <div style="display: flex">
    <textarea rows="1" wrap="off" readonly style="display: inline-block; width: 100%; resize: none; overflow-y: hidden"><a href="https://www.rdkit.org/"><img alt="Powered by RDKit" src="https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC"></a>
    </textarea>
    <button onclick="copyBadgeMarkdown(this)" style="display: inline-block; padding: 0 8px">Copy</button>
  </div>
</details>

<script>
  function copyBadgeMarkdown(btn) {
    if (navigator.clipboard) {
      var text = btn.previousElementSibling.textContent;
      navigator.clipboard.writeText(text);
    }
  }
</script>

## Integration with other open-source projects
- [KNIME](https://www.knime.com/rdkit): Workflow and analytics tool
- [PostgreSQL](https://github.com/rdkit/rdkit/blob/master/Docs/Book/Cartridge.rst): Extensible relational database
- [Django](http://django-rdkit.readthedocs.org/en/latest/): "The web framework for perfectionists with deadlines"
- [SQLite](https://github.com/rvianello/chemicalite): "The most used database engine in the world"
- [Lucene](https://github.com/rdkit/org.rdkit.lucene): Text-search engine [^footnote1]

## The Contrib Directory

The Contrib directory, part of the standard RDKit distribution, includes code that has been contributed by members of the community.


[^footnote1]: These implementations are functional but are not necessarily the best, fastest, or most complete.


## License

This document is copyright (C) 2013-2024 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit <http://creativecommons.org/licenses/by-sa/4.0/> or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
