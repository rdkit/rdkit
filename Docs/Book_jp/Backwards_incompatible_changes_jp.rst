後方互換性の無い変更
################################
[`Backwards incompatible changes <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#backwards-incompatible-changes>`__]

リリース間の後方互換性を維持しようと本当に努めていますが、互換性を破る必要がある事例がまれにあります。
このドキュメントの目的は、RDKitで行われた後方互換性の無い変更についての情報を補足することです。
リリースサイクルによって分類されており、だいたいにおいてバグの修正で生じた結果による変更は\ *含んでいません*\ 。
バグの修正はリリースノートで注意喚起するようにしています。

Release 2018.03
**********************************
[`Release 2018.03 <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#release-2018-03>`__]

``MolToSmiles()``\ がデフォルトでisomeric SMILESを生成
===========================================================================================
[`MolToSmiles() generates isomeric SMILES by default <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#moltosmiles-generates-isomeric-smiles-by-default>`__]

以前のリリースでは、出力されるSMILESに立体化学や同位体ラベルについての情報を含めたい場合、
オプションの\ ``isomericSmiles``\ 引数をtrueに設定する必要がありました。
現在、このデフォルトの値はtrueとなっています。昔の動作に戻し、立体化学の情報を持たないSMILESを取得したい場合は、\ ``isomericSmiles``\ をfalseに設定するだけです。

``MolToMolBlock()``\ はincludeStereoフラグがセットされると2Dコンフォメーションを生成する
===========================================================================================
[`MolToMolBlock() generates a 2D conformation when the includeStereo flag is set <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#moltomolblock-generates-a-2d-conformation-when-the-includestereo-flag-is-set>`__]

Mol blockの立体化学を取得したい場合、出力が座標の情報を有する必要があります。
以前のRDKitのバージョンでは、座標の生成を忘れずに自分で行う必要がありました。現在では、\ ``includeStereo``\ がセットされている場合、コンフォメーションを持たない分子に対してはデフォルトで座標が生成されます。

現在のデフォルトではコンフォメーション生成コードはETKDGを使う
===========================================================================================

[`The conformation generation code now uses ETKDG by default <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#the-conformation-generation-code-now-uses-etkdg-by-default>`__]

以前のRDKitのリリースでは、デフォルトでは、標準的なデジスタンスジオメトリー法を使ってコンフォメーションを生成していました。新しいデフォルトの設定では、Sereina
RinikerによるETKDGアルゴリズムを使います。少し遅いですが、ずっと良い結果を出すことが示されています。

Release 2018.09
***********************************
[`Release 2018.09 <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#release-2018-09>`__]

``GetAtomSmile()``\ はデフォルトでisomeric SMILESを生成
===========================================================================================
[`GetAtomSmiles() generates isomeric SMILES by default <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#getatomsmiles-generates-isomeric-smiles-by-default>`__]

以前のリリースでは、出力されるSMILESに立体化学や同位体ラベルについての情報を含めたい場合、オプションの\ ``isomericSmiles``\ 引数をtrueに設定する必要がありました。
現在、このデフォルトの値はtrueとなっています。昔の動作に戻し、立体化学の情報を持たないSMILESを取得したい場合は、\ ``isomericSmiles``\ をfalseに設定するだけです。

MCSコードのringMatchesRingOnlyオプションの変更点
===========================================================================================
[`Changes to ringMatchesRingOnly option in the MCS code <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#changes-to-ringmatchesringonly-option-in-the-mcs-code>`__]

FindMCS()関数のringMatchesRingOnlyオプションは結合と結合のマッチングと同様に、原子と原子のマッチングにも適用されます。

現在のデフォルトではPythonから呼び出すとコンフォメーション生成コードはETKDGを使う
===========================================================================================
[`The conformation generation code now uses ETKDG by default when called from Python <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#the-conformation-generation-code-now-uses-etkdg-by-default-when-called-from-python>`__]

現在デフォルトでは、Python関数のEmbedMolecule()とEmbedMultipleConfs()は標準的なデジスタンスジオメトリーの代わりにETKDGアルゴリズムを使用します。

ライセンス
*****************
[`License <https://www.rdkit.org/docs/BackwardsIncompatibleChanges.html#license>`__]

この文書の著作権は copyright (C) 2013-2018 by Greg Landrum
に所属しています。

この文書はCreative Commons Attribution-ShareAlike 4.0
Licenseのもとでライセンスされています。このライセンスを見るためには http://creativecommons.org/licenses/by-sa/4.0/
にアクセスするか、Creative Commons, 543 Howard Street, 5th Floor, San
Francisco, California, 94105, USA.　に手紙を送ってください。

このライセンスの意図はRDKitそのものの意図と似ています。簡単に言えば
“これを使ってなんでもやりたいことをやっていいですが、私たちの功績についても言及してください”
