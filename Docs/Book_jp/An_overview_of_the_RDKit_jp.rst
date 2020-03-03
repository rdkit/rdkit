RDKitの概要
#################################################################
[`An overview of the RDKit <https://www.rdkit.org/docs/Overview.html#an-overview-of-the-rdkit>`__]

RDKitっていったい？
***************************************************************
[`What is it? <https://www.rdkit.org/docs/Overview.html#what-is-it>`__]

ケモインフォマティクスのためのオープンソースのツールキット
===========================================================================
[`Open source toolkit for cheminformatics <https://www.rdkit.org/docs/Overview.html#open-source-toolkit-for-cheminformatics>`__]

-  ビジネスフレンドリーなBSDライセンス
-  コアデータストラクチャーとアルゴリズムはC++で記載
-  Python(2.x および3.x)のラッパーは Boost.Pythonを使って作成
-  Java と C# のラッパーはSWIGを使って作成
-  2D と 3D で分子を取り扱えます
-  機械学習のための記述子を生成
-  PostgreSQLのための化合物データベース・カートリッジ
-  KNIMEのためのケモインフォマティクスノード (KNIMEコミュニティサイト：https://www.knime.com/rdkit から配布されています)

運営に関する情報
==========================================================================
[`Operational <https://www.rdkit.org/docs/Overview.html#operational>`__]

-  http://www.rdkit.org
-  サポート　Mac/Windows/Linux 　
-  6ヶ月ごとにリリース
-  Web上の情報：

   -  ホームページ： http://www.rdkit.org ドキュメンテーション、リンク
   -  Github:(https://github.com/rdkit)　ダウンロード、バグ追跡、Gitリポジトリ
   -  Sourceforge (http://sourceforge.net/projects/rdkit) メーリングリスト
   -  ブログ (https://rdkit.blogspot.com) うまく使うための秘訣, 技, 雑多な話題
   -  チュートリアル (https://github.com/rdkit/rdkit-tutorials)
      RDKit使用方法のJupyterベースのチュートリアル
   -  KNIMEとの統合 (https://github.com/rdkit/knime-rdkit)
      KNIME向けのRDKit ノード
   -  メーリングリストは https://sourceforge.net/p/rdkit/mailman/
      で、検索可能なアーカイブは\ `rdkit-discuss <https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/>`__\ と\ `rdkit-devel <https://www.mail-archive.com/rdkit-devel@lists.sourceforge.net/>`__\ です

-  ソーシャルメディア：

   -  Twitter: @RDKit_org
   -  LinkedIn: https://www.linkedin.com/groups/8192558
   -  Slack: https://rdkit.slack.com
      (招待が必要です、Gregに連絡してください）

沿革
=================================================
[`History <https://www.rdkit.org/docs/Overview.html#history>`__]

-  2000-2006: Rational
   DiscoveryにおいてADME、Tox、生理活性を予測するモデルを作成するために開発され使われていた
-  June 2006: ソフトウェアをオープンソース(BSD
   license)としてリリース、Rational Discoveryを閉鎖
-  現在まで:
   オープンソースでの開発を継続、Novartis内で使用、Novartisからのコントリビューションはオープンソース版へ反映

他のオープンソースプロジェクトとの統合
****************************************************************
[`Integration with other open-source projects <https://www.rdkit.org/docs/Overview.html#integration-with-other-open-source-projects>`__]

-  `KNIME <https://www.knime.com/rdkit>`__:ワークフローと分析ツール
-  `PostgreSQL <https://www.rdkit.org/docs/Cartridge.html>`__: 拡張可能なリレーショナルデータベース
-  `Django <http://django-rdkit.readthedocs.org/en/latest/>`__
   :“締め切りに追われる完璧主義者のためのウェブフレームワーク”
-  `SQLite <https://github.com/rvianello/chemicalite>`__
   -“世界で最も使われているデータベースエンジン”
-  `Lucene <https://github.com/rdkit/org.rdkit.lucene>`__ :
   テキストサーチエンジン  [#f1]_

.. 訳注 PostgreSQLリンク先が https://github.com/rdkit/rdkit/blob/master/Docs/Book/Cartridge.rst となっており機能しなかったためRDKit database cartridgeのリンクとした。

他のオープソースプロジェクトによる使用
*****************************************************************************
[`Usage by other open-source projects <https://www.rdkit.org/docs/Overview.html#usage-by-other-open-source-projects>`__]

この項目は必然的に最新の情報ではなくなってしまいます。もし他にご存知のプロジェクトがあれば我々に連絡、あるいはプルリクエストをサブミットしてください！

-  `gpusimilarity <https://github.com/schrodinger/gpusimilarity>`__ -
   フィンガープリント類似性探索のCuda/Thrustによる実装
-  `Samson Connect <https://www.samson-connect.net/>`__ -
   ナノシステムのシミュレーションと適合モデリングのためのソフトウェア
-  `mol_frame <https://github.com/apahl/mol_frame>`__ - DaskとPandas
   DataFrameのための化学構造の取り扱い 　
-  `RDKitjs <https://github.com/cheminfo/RDKitjs>`__ -
   JavasScriptのためのRDKit機能のポート
-  `DeepChem <https://deepchem.io/>`__ -
   化学のためのディープラーニングのためのpythonライブラリ
-  `mmpdb <https://github.com/rdkit/mmpdb>`__ - Matched molecular pair
   データベースの生成と分析
-  `CheTo <https://github.com/rdkit/CheTo>`__
   (`文献 <http://pubs.acs.org/doi/10.1021/acs.jcim.7b00249>`__)-
   ケミカルトピックモデリング
-  `OCEAN <https://github.com/rdkit/OCEAN>`__\ (`文献 <http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00067>`__)-
   最適化交差反応性予測
-  `ChEMBL Beaker <https://github.com/mnowotka/chembl_beaker>`__ - RDKitとOSARのためのスタンドアローンウェブサーバーラッパー
-  `myChEMBL <https://github.com/chembl/mychembl>`__
   (`ブログ記事 <http://chembl.blogspot.de/2013/10/chembl-virtual-machine-aka-mychembl.html>`__,
   `文献 <http://bioinformatics.oxfordjournals.org/content/early/2013/11/20/bioinformatics.btt666>`__
   ) -ケモインフォマティクスツールとオープンデータの仮想マシンによる実装
-  `ZINC <http://zinc15.docking.org/>`__ -
   バーチャルスクリーニングのための購入可能な化合物の無料データベース
-  `sdf_viewer.py <https://github.com/apahl/sdf_viewer>`__ -
   インラタクティブなSDFビューワー
-  `sdf2ppt <https://github.com/dkuhn/sdf2ppt>`__ -
   Powerpoint/openofficeプレゼンテーションでSDFファイルの読み込みとイメージグリッドとして分子の表示を行う
-  `MolGears <https://github.com/admed/molgears>`__ -
   生理活性分子のためのケモインフォマティクスツール
-  `PYPL <http://www.biochemfusion.com/downloads/#OracleUtilities>`__ -
   Oracle PL/SQLからPythonスクリプトの呼び出しを可能にするシンプルなカートリッジ
-  `shape-it-rdkit <https://github.com/jandom/shape-it-rdkit>`__ -
   silicos itからRDKitバックエンドに移植されたガウス型の分子重ね合わせコード
-  `WONKA <http://wonka.sgc.ox.ac.uk/WONKA/>`__ -
   タンパク質-リガンド共結晶構造の解析と精査のためのツール
-  `OOMMPPAA <http://oommppaa.sgc.ox.ac.uk/OOMMPPAA/>`__ -
   タンパク質-リガンド共結晶構造に基づく指向性合成とデータ分析のためのツール
-  `OCEAN <https://github.com/rdkit/OCEAN>`__ -
   ChEMBLをデータソースとして用いて化学構造の標的を予測するためのツール
-  `chemfp <http://chemfp.com/>`__ - 非常に高速なフィンガープリント検索
-  `rdkit_ipynb_tools <https://github.com/apahl/rdkit_ipynb_tools>`__ -
   IPython NotebookのためのRDKitツール
-  `Vernalis KNIME nodes <https://www.knime.com/book/vernalis-nodes-for-knime-trusted-extension>`__
-  `Erlwood KNIME nodes <https://www.knime.com/community/erlwood>`__
-  `AZOrange <https://github.com/AZcompTox/AZOrange>`__

Contribディレクトリ
******************************************************************
[`The Contrib Directory <https://www.rdkit.org/docs/Overview.html#the-contrib-directory>`__]

標準RDKitディストリビューションの一部である Contrib
ディレクトリはコミュニティメンバーからのコントリビューションによるコードを含みます。

脚注
*****************************************************************
[`Footenotes <https://www.rdkit.org/docs/Overview.html#footnotes>`__]

.. rubric:: Footnotes

.. [#f1] これらの実装は機能しますが、必ずしも最良、最速、あるいは最も完全であるとは限りません。

ライセンス
**************************************************************
[`License <https://www.rdkit.org/docs/Overview.html#license>`__]

この文書の著作権は copyright (C) 2013-2018 by Greg Landrum
に所属しています。

この文書はCreative Commons Attribution-ShareAlike 4.0
Licenseのもとでライセンスされています。このライセンスを見るためには http://creativecommons.org/licenses/by-sa/4.0/
にアクセスするか、Creative Commons, 543 Howard Street, 5th Floor, San
Francisco, California, 94105, USA.　に手紙を送ってください。

このライセンスの意図はRDKitそのものの意図と似ています。簡単に言えば
“これを使ってなんでもやりたいことをやっていいですが、私たちの功績についても言及してください”
