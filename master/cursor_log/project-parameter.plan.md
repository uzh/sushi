# プロジェクトごとのデフォルトパラメータ機能の実装

## 概要

各プロジェクトフォルダ（例：`/srv/gstore/projects/p1234/`）にデフォルトパラメータファイル（`project_default_parametersets.tsv`）を保存し、ジョブ実行時にパラメータを以下の優先順位で適用する機能を実装します：

1. 親データセットから継承されたパラメータ（例：strandMode, refBuild）
2. プロジェクトデフォルトファイルから読み込まれたパラメータ
3. アプリケーションのデフォルトパラメータ

パラメータはアプリケーション名でスコープ化され（例：`FastqcApp::cores`）、異なるアプリ間での競合を防ぎます。

## 主要な変更点

### 1. SushiAppクラスの拡張（`lib/sushi_fabric/lib/sushi_fabric/sushiApp.rb`）

- **プロジェクトデフォルトファイルの読み込みメソッド** (`load_project_defaults`)
  - プロジェクトルートから`project_default_parametersets.tsv`を読み込む
  - App名でスコープ化されたパラメータを解析（例：`FastqcApp::cores`）
  - 現在のアプリに関連するパラメータのみを抽出

- **プロジェクトデフォルトファイルの保存メソッド** (`save_project_defaults`)
  - 既存のファイルを読み込んで既存パラメータを保持
  - 現在のアプリのパラメータを更新または追加
  - App名でスコープ化してファイルに書き込む

- **`set_default_parameters`メソッドの拡張**
  - 現在のロジック：アプリのデフォルトを設定
  - 追加：プロジェクトデフォルトを読み込んで適用（親データセットパラメータを上書きしない）

### 2. データセット画面の更新（`app/views/data_set/_sushi_application_list.html.erb`）

- アプリケーションボタンフォームに`use_project_defaults`チェックボックスを追加
  - 配置：「Applications - refresh」の下、アプリケーションボタンテーブルの前
  - ラベル：「Use project default parameters」
  - デフォルト：checked（プロジェクトデフォルトを使用する）
  - フォームの一部としてパラメータを`set_parameters`アクションに渡す

### 3. 確認画面の更新（`app/views/run_application/confirmation.html.erb`）

- `save_as_default`チェックボックスを追加
  - ラベル：「Save these parameters as project defaults」
  - デフォルト：unchecked
  - Submit/MockRunボタンの近くに配置
  - 目的：ジョブサブミット時にプロジェクトデフォルトファイルに保存

### 4. コントローラーの更新（`app/controllers/run_application_controller.rb`）

- **`set_parameters`メソッドを更新**
  - `params[:use_project_defaults]`をチェック
  - trueの場合、`@sushi_app.set_default_parameters`の後に`@sushi_app.load_project_defaults`を呼び出す
  - これによりプロジェクトデフォルトが既存パラメータに上書きされる

- **`submit_jobs`メソッドを更新**
  - `save_as_default`パラメータを`active_job_params`に追加
  - パラメータをジョブに渡す

### 5. ジョブ実行時の保存処理（`app/jobs/submit_job.rb`）

- `perform`メソッドを更新
  - `save_as_default`フラグをチェック
  - trueの場合、`sushi_app.save_project_defaults`を呼び出してパラメータを保存
  - ジョブ実行前または実行後に保存（エラー時にも保存されるように実行前が望ましい）

## ファイル形式

`project_default_parametersets.tsv`の形式：

```
FastqcApp::cores	8
FastqcApp::ram	16
STARApp::refBuild	Homo_sapiens/Ensembl/GRCh38
STARApp::cores	16
```

- タブ区切り形式
- 1行目：`AppName::parameter_name<TAB>value`
- スコープ化により異なるアプリで同じパラメータ名があっても競合しない

## パラメータ優先順位の実装

`set_default_parameters`メソッドで以下の順序で適用：

1. アプリのデフォルトパラメータを設定（各Appのinitializeとset_default_parameters）
2. プロジェクトデフォルトを読み込んで適用（既存値がない場合のみ）
3. 親データセットのパラメータを適用（最優先、常に上書き）

コントローラー側（`set_parameters`）では、親データセットのパラメータを明示的に適用する処理を追加する必要がある場合があります。

## 注意点

- プロジェクトフォルダが存在しない場合はエラーを出さずスキップ
- ファイルが存在しない場合も静かに処理を続行
- ファイルの読み書き時に適切なエラーハンドリング
- プロジェクトパス（`@project`）は既に`p1234`形式で持っているので、`@gstore_project_dir`を使用

## To-dos

- [ ] SushiAppクラスにload_project_defaultsメソッドを追加してプロジェクトデフォルトファイルを読み込む
- [ ] SushiAppクラスにsave_project_defaultsメソッドを追加してプロジェクトデフォルトファイルを保存する
- [ ] set_default_parametersの後にload_project_defaultsを呼び出してプロジェクトデフォルトを適用
- [ ] データセット画面にuse_project_defaultsチェックボックスを追加
- [ ] run_application_controllerのset_parametersでuse_project_defaultsをチェックして適用
- [ ] confirmation.html.erbにsave_as_defaultチェックボックスを追加
- [ ] run_application_controllerのsubmit_jobsでsave_as_defaultパラメータをactive_job_paramsに渡す
- [ ] submit_job.rbでsave_as_defaultがtrueの場合にプロジェクトデフォルトを保存

