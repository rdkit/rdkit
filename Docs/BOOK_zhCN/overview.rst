.. _overview:

==================
RDKit简介
==================
RDKt 是非常强大的化学信息python工具包。
[`An overview of the RDKit <https://www.rdkit.org/docs/Overview.html#an-overview-of-the-rdkit>`__]


在2019-ncov期间，2020-02-03 对该文档进行了翻译！ 


RDKit是什么？
***************************************************************

开源的化学信息工具包
===========================================================================
-  采用了商业友好的BSD协议
-  核心数据结构和算法用C++实现
-  通过Boost.Python技术对RDKit进行封装，提供Python2/Python3的接口 
-  通过SWIG技术提供Java 和C# 接口
-  提供了大量对化学分子2D/3D的计算操作
-  生成用于机器学习的分子描述符
-  基于PostgreSQL搭建分子数据库 
-  KNIME中的化学信息计算支持（https://www.knime.com/rdkit）

RDKit其他信息
==========================================================================

- RDKit官网： http://www.rdkit.org
- 支持多种系统：　Mac/Windows/Linux 　
-  版本发布周期： 6个月 
-  Web上的资源

   -  官网： http://www.rdkit.org 
   -  Github:(https://github.com/rdkit)　
   -  Sourceforge (http://sourceforge.net/projects/rdkit) 
   -  博客： (https://rdkit.blogspot.com)  
   -  教程： (https://github.com/rdkit/rdkit-tutorials)
   -  KNIME和RDKit的集成 (https://github.com/rdkit/knime-rdkit)
- 邮件列表： https://sourceforge.net/p/rdkit/mailman/
-  社交媒体：

   -  Twitter: @RDKit_org
   -  LinkedIn: https://www.linkedin.com/groups/8192558
   -  Slack: https://rdkit.slack.com（需要加入的话，联系Greg）
- 教程：
   - 1. “Getting Started width the RDKit in Python”, 
   - Greg Landrum’s RDkit blog,
   - and a particularly nice example from the ChEMBL-og. 

发展过程
=================================================

-  2000-2006: 为了建立ADMET、生物活性等性质的预测模型，开发了相应的工具。
-  June 2006: 发布该工具（BSD license）
-  至今: 诺华内部支持继续更新和完善该软件，并对社区开源。

和其他开源工具整合
****************************************************************

-  `KNIME <https://www.knime.com/rdkit>`__:工作流程和分析工具
-  `PostgreSQL <https://www.rdkit.org/docs/Cartridge.html>`__: 可拓展关系数据库
-  `Django <http://django-rdkit.readthedocs.org/en/latest/>`__：python的web框架
   :“締め切りに追われる完璧主義者のためのウェブフレームワーク”
-  `SQLite <https://github.com/rvianello/chemicalite>`__：使用最多的数据库，最简单的数据库
-  `Lucene <https://github.com/rdkit/org.rdkit.lucene>`__ :文本搜索引擎[#f1]_


使用RDKit的开源工具 
*****************************************************************************
- `Open Force Field Toolkit <https://github.com/openforcefield/openforcefield/>'__ - 于直接化学感知的力场参数引擎
- `stk <https://github.com/lukasturcani/stk>`__  - 自动设计分子的python包
-  `gpusimilarity <https://github.com/schrodinger/gpusimilarity>`__ - 基于gpu的分子指纹搜索工具
-  `Samson Connect <https://www.samson-connect.net/>`__ - 纳米系统建模仿真软件
-  `mol_frame <https://github.com/apahl/mol_frame>`__ - 为Dask和Pandas 提供处理化学结构的接口 
-  `RDKitjs <https://github.com/cheminfo/RDKitjs>`__ -  JavasScript封装了RDKit的部分功能
-  `DeepChem <https://deepchem.io/>`__ - 用于化学领域深度学习的python包
-  `mmpdb <https://github.com/rdkit/mmpdb>`__ - Matched molecular pair database generation and analysis 
-  `CheTo <https://github.com/rdkit/CheTo>`__ - Chemical topic modeling
-  `OCEAN <https://github.com/rdkit/OCEAN>`__ - Optimized cross reactivity estimation
-  `ChEMBL Beaker <https://github.com/mnowotka/chembl_beaker>`__ - RDKit 和 OSRA 的独立 web 服务 
-  `myChEMBL <https://github.com/chembl/mychembl>`__ - A virtual machine implementation of open data and cheminformatics tools
-  `ZINC <http://zinc15.docking.org/>`__ - 用于虚拟筛选的可购买化合物开源数据库
-  `sdf_viewer.py <https://github.com/apahl/sdf_viewer>`__ - sdf 查看器
-  `sdf2ppt <https://github.com/dkuhn/sdf2ppt>`__ - 把sdf转换成PPT
-  `MolGears <https://github.com/admed/molgears>`__ - A cheminformatics tool for bioactive molecules
-  `PYPL <http://www.biochemfusion.com/downloads/#OracleUtilities>`__ - Simple cartridge that lets you call Python scripts from Oracle PL/SQL.
-  `shape-it-rdkit <https://github.com/jandom/shape-it-rdkit>`__ -  Gaussian molecular overlap code shape-it (from silicos it) ported to RDKit backend
-  `WONKA <http://wonka.sgc.ox.ac.uk/WONKA/>`__ -  蛋白配体复合物结构分析和查询的服务
-  `OOMMPPAA <http://oommppaa.sgc.ox.ac.uk/OOMMPPAA/>`__ - Tool for directed synthesis and data analysis based on protein-ligand crystal structures
-  `OCEAN <https://github.com/rdkit/OCEAN>`__ - web-tool for target-prediction of chemical structures which uses ChEMBL as datasource
-  `chemfp <http://chemfp.com/>`__ - 快速分子指纹相似性搜索 
-  `rdkit_ipynb_tools <https://github.com/apahl/rdkit_ipynb_tools>`__ - 用于jupyter notebook的rdkit工具
-  `Vernalis KNIME nodes <https://www.knime.com/book/vernalis-nodes-for-knime-trusted-extension>`__
-  `Erlwood KNIME nodes <https://www.knime.com/community/erlwood>`__
-  `AZOrange <https://github.com/AZcompTox/AZOrange>`__

Contrib目录
******************************************************************
contrib目录位置： https://github.com/rdkit/rdkit/tree/master/Contrib 。
contrib目录，是RDKit发布版本的一部分，包含社区成员贡献的代码。

脚注
*****************************************************************
.. rubric:: Footnotes

.. [#f1] These implementations are functional but are not necessarily the best, fastest, or most complete.

协议
**************************************************************

This document is copyright (C) 2013-2018 by Greg Landrum

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License. To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.

The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”



