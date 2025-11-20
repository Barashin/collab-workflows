# ワークフロー各ノードの説明

このディレクトリには、ワークフロー画像の各ノードに対応するPythonスクリプトが含まれています。

## ノード一覧

### Node 1: リガンドファイルダウンロード
- **ファイル**: `ligand/node_01_download_ligands.py`
- **入力**: なし（暗黙的）
- **出力**: `ligands.zip`
- **説明**: ZenodoからリガンドライブラリのZIPファイルをダウンロード

### Node 2: PDBファイルダウンロード
- **ファイル**: `protein/node_05_download_pdb.py`
- **入力**: `PDB_ID` (環境変数またはデフォルト値 "5Y7J")
- **出力**: `{pdb_id}.pdb` (例: `5Y7J.pdb`)
- **説明**: RCSB PDBからタンパク質構造ファイルをダウンロード

### Node 3: リガンドファイル展開
- **ファイル**: `ligand/node_02_unpack_ligands.py`
- **入力**: `ligands.zip`
- **出力**: `ligand/output/constructed_library/` ディレクトリ（複数のSDFファイル）
- **説明**: ZIPファイルを展開してリガンドSDFファイルを取得

### Node 4: タンパク質入力確認
- **ファイル**: `protein/node_06_protein_input.py`
- **入力**: `PDB_ID` (環境変数またはデフォルト値)
- **出力**: `{pdb_id}.pdb` (確認のみ)
- **説明**: PDBファイルが存在することを確認

### Node 5: リガンド選択
- **ファイル**: `ligand/node_03_ligand_selection.py`
- **入力**: `ligand/output/constructed_library/` ディレクトリ内のSDFファイル群
- **出力**: `ligand/output/Ligands_select.sdf`
- **説明**: 選択されたリガンドを1つのSDFファイルにまとめる
- **注意**: デフォルトではテスト用に `clean_drug108*.sdf` パターンのファイルを選択

### Node 6: リガンドビュー
- **ファイル**: `ligand/node_10_ligand_view.py`
- **入力**: `ligand/output/Ligands_select.sdf`
- **出力**: `ligand/output/ligand.csv`
- **説明**: リガンドの情報（SMILES、分子量、LogPなど）をCSV形式で出力

### Node 7: Biopython - チェーン抽出
- **ファイル**: `preparation/node_08_extract_chains.py`
- **入力**: `{pdb_id}.pdb`
- **出力**: `preparation/output/{pdb_id}_chain.pdb`
- **説明**: PDBファイルからA鎖とB鎖、およびリファレンスリガンドを抽出

### Node 8: Biopython - リガンド中心識別
- **ファイル**: `preparation/node_11_ligand_center_identification.py`
- **入力**: `ligand/output/Ligands_select.sdf`, `ligand/output/ligand.csv`, `preparation/output/{pdb_id}_chain.pdb`
- **出力**: `preparation/output/config.txt`
- **説明**: リガンドの中心座標を計算し、ドッキング設定ファイルを生成

### Node 9: OpenMM PDBFixer - タンパク質クリーンアップ
- **ファイル**: `preparation/node_12_protein_extraction.py`
- **入力**: `preparation/output/{pdb_id}_chain.pdb`
- **出力**: `preparation/output/{pdb_id}_clean.pdb`
- **説明**: PDBFixerを使用してタンパク質構造をクリーンアップ（欠損残基の追加、水素の追加など）

### Node 10: PDB2PQR - AMBER電荷適用
- **ファイル**: `preparation/node_10_apply_charges.py` (存在する場合)
- **入力**: `preparation/output/{pdb_id}_clean.pdb`
- **出力**: `preparation/output/{pdb_id}_amber.pqr`, `preparation/output/{pdb_id}_amber.pdb`
- **説明**: PDB2PQRを使用してAMBER力場の電荷を適用

### Node 11: SMINA - インシリコスクリーニング
- **ファイル**: `docking/node_13_smina_screening.py`
- **入力**: 
  - `preparation/output/config.txt`
  - `preparation/output/{pdb_id}_clean.pdb` (または `{pdb_id}_amber.pqr`)
  - `ligand/output/selected_compounds/` (個別のSDFファイル) または `ligand/output/Ligands_select.sdf`
- **出力**: `docking/output/` ディレクトリ（ドッキング結果ファイル群）
- **説明**: SMINAを使用して分子ドッキングを実行

### Node 12: レポート生成
- **ファイル**: `report/node_14_reporting.py` または `report.py` (統合スクリプト)
- **入力**: `docking/output/` ディレクトリ内のドッキング結果
- **出力**: 
  - `report/output/docking_ranking.txt` または `results/docking_ranking.txt`
  - `report/output/*_docked.sdf` (トップ化合物)
- **説明**: ドッキング結果を解析してランキングを生成

## 統合スクリプト

### `protein_preparation.py`
複数のノード処理を統合したスクリプト：
- PDBファイルのダウンロード
- AB鎖の抽出
- リガンドの中心座標計算
- グリッドサイズ計算とconfigファイル生成
- PDBFixerによる構造修正
- PDB2PQRによる電荷付与

### `report.py`
レポート生成の統合スクリプト：
- ドッキング結果の解析とランキング生成
- トップ化合物のファイルコピー

## 実行順序

ワークフローは以下の順序で実行する必要があります：

```
1. Node 1 (リガンドダウンロード) ─┐
2. Node 2 (PDBダウンロード)      ─┤ 並列実行可能
                                   │
3. Node 3 (リガンド展開) ← Node 1
4. Node 4 (タンパク質入力確認) ← Node 2
5. Node 5 (リガンド選択) ← Node 3
6. Node 6 (リガンドビュー) ← Node 5
7. Node 7 (チェーン抽出) ← Node 2
                                   │
8. Node 8 (リガンド中心識別) ← Node 5, 6, 7
9. Node 9 (タンパク質クリーンアップ) ← Node 7
10. Node 10 (AMBER電荷適用) ← Node 9
                                  │
11. Node 11 (SMINAスクリーニング) ← Node 8, 10
12. Node 12 (レポート生成) ← Node 11
```

## 実行方法

### 一括実行（推奨）

ワークフロー全体を3つのステージに分けて実行：

```bash
# ステージ1: データ準備（タンパク質のダウンロードと準備）
bash pre_run.sh

# ステージ2: ドッキングスクリーニング
bash run.sh

# ステージ3: レポート生成
bash post_run.sh
```

### 統合スクリプト実行

統合されたスクリプトを使用する場合：

```bash
# タンパク質準備（PDBダウンロード、チェーン抽出、構造修正、電荷付与を含む）
python protein_preparation.py

# ドッキングスクリーニング
python docking/node_13_smina_screening.py

# レポート生成
python report.py
```

### 個別実行

各ノードを個別に実行する場合：

```bash
# リガンド関連
python ligand/node_01_download_ligands.py
python ligand/node_02_unpack_ligands.py
python ligand/node_03_ligand_selection.py
python ligand/node_04_prepare_ligands.py
python ligand/node_09_real_ligand_addition.py
python ligand/node_10_ligand_view.py

# タンパク質関連
python protein/node_05_download_pdb.py
python protein/node_06_protein_input.py

# 準備関連
python preparation/node_07_ligand_loading.py
python preparation/node_08_extract_chains.py
python preparation/node_11_ligand_center_identification.py
python preparation/node_12_protein_extraction.py

# ドッキング
python docking/node_13_smina_screening.py

# レポート
python report/node_14_reporting.py
# または
python report.py
```

## 環境変数

- `PDB_ID`: PDB IDを指定（デフォルト: "5Y7J"）

```bash
export PDB_ID="5Y7J"
```

## ファイル構造

```
updated_job-1/
├── pre_run.sh              # データ準備ステージ
├── run.sh                  # ドッキングスクリーニングステージ
├── post_run.sh             # レポート生成ステージ
├── protein_preparation.py  # 統合タンパク質準備スクリプト
├── report.py               # 統合レポート生成スクリプト
├── ligand/                 # リガンド処理ノード
│   ├── node_01_download_ligands.py
│   ├── node_02_unpack_ligands.py
│   ├── node_03_ligand_selection.py
│   ├── node_04_prepare_ligands.py
│   ├── node_09_real_ligand_addition.py
│   ├── node_10_ligand_view.py
│   └── output/             # リガンド処理の出力
├── protein/                # タンパク質処理ノード
│   ├── node_05_download_pdb.py
│   └── node_06_protein_input.py
├── preparation/            # 準備処理ノード
│   ├── node_07_ligand_loading.py
│   ├── node_08_extract_chains.py
│   ├── node_11_ligand_center_identification.py
│   ├── node_12_protein_extraction.py
│   └── output/             # 準備処理の出力
├── docking/                # ドッキングノード
│   ├── node_13_smina_screening.py
│   └── output/             # ドッキング結果
└── report/                 # レポートノード
    ├── node_14_reporting.py
    └── output/             # レポート出力
```

## 注意事項

1. **Node 5 (リガンド選択)**: デフォルトではテスト用に `clean_drug108*.sdf` パターンのファイルを選択します。全ライブラリを処理する場合は、`ligand/node_03_ligand_selection.py` の `LIGAND_PATTERN` を変更してください。

2. **Node 11 (SMINAスクリーニング)**: `ligand/output/Ligands_select.sdf` が存在しない場合、`ligand/output/selected_compounds/` ディレクトリ内の個別のSDFファイルを直接処理します。

3. **依存関係**: 各ノードは前のノードの出力を入力として使用するため、順序を守って実行してください。

4. **統合スクリプト**: `protein_preparation.py` は複数のノード処理を統合したスクリプトです。個別のノードを実行する代わりに、このスクリプトを使用することもできます。

5. **出力ディレクトリ**: 
   - `ligand/output/`: リガンド処理の出力
   - `preparation/output/`: タンパク質準備の出力
   - `docking/output/`: ドッキング結果
   - `report/output/` または `results/`: 最終結果（ランキングと最良のリガンド）
