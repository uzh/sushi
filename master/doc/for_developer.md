# SUSHI System Developer Document

## 1. System Overview
### 1.1 What is SUSHI System
- SUSHI System is a web application for genome data analysis, developed based on Rails.
- Users can upload data, execute analysis jobs, and retrieve results.
- The Analysis Pipeline is managed by **SUSHIApp** and executed through **ezRun**.
- **SUSHIApp corresponds to a single analysis pipeline and generates both analysis results and metadata. The metadata associated with SUSHIApp is called DataSet, which serves as both input and output for SUSHIApp (it is actually a flexible two-dimensional table serialized as a Hash object and stored as a text field in the database).**

### 1.2 System Architecture
- **Backend**: Ruby on Rails
- **Frontend**: ERB / JavaScript (Uses JQuery. Vue.js and React are not used.)
- **Data Management**: MySQL (Manages metadata and all non-analysis result data), gStore (Stores analysis data at /srv/gstore/)
- **Analysis Pipeline**: SUSHIApp (Stored as `.rb` files under `lib/` directory), ezRun (R library)
- **Authentication & Authorization**: Devise + Center-internal LDAP authentication
- **Job Management**: Job Manager (Module that bridges to Slurm)
- **Cache Management**: KeyDB (Used for dataset tree caching)

---

## 2. Core System Components
### 2.1 Data Management
- **gStore**: Storage for analysis data (/srv/gstore/)
- **MySQL**: Manages analysis metadata and all other non-analysis result data

### 2.2 Analysis Pipeline
- **SUSHIApp**
  - Each analysis pipeline is stored as a `.rb` file under the `lib/` directory.
  - SUSHIApp calls **ezRun** to execute R scripts.

### 2.3 Authentication & Authorization
- **Devise** is used for user management.
- **Authentication is delegated to the center-internal LDAP server.**

### 2.4 Job Management
- **Job Manager**
  - The use of `job_manager.rb` in Workflow Manager and Sidekiq is currently discontinued.
  - Acts as a bridge to **Slurm** for managing analysis jobs.
  - Jobs are submitted via Rails.
- **Workflow Manager**
  - Uses **KeyDB** for dataset tree caching.

---

## 3. Development Workflow
### 3.1 Test SUSHI Setup
- Set up the SUSHI development environment and verify the operation of each component.
- SUSHI の開発環境をセットアップし、各コンポーネントの動作を確認する。

### 3.2 Test Run Procedure
- Apply development changes and execute a test job.
- Check the logs to ensure expected operation.


### 3.3 Commit Procedure
- Commit changes to Git following the appropriate branching strategy.
- Create a Pull Request and merge after code review.

---

## 4. Deployment & Maintenance
### 4.1 Deployment Procedure (ToDo)
- **This section needs to be reviewed and revised.**。

---

## 5. Known Issues & Future Development
### 5.1 Known Issues and Limitations
- **No unit tests available.**
- **No CI/CD.** The custom setup script (shell script) only handles SUSHI server type-specific code modifications.
- **UI is based on legacy code from the Rails 4 era, and Hotwire (Turbo & Stimulus) is not supported.** → **A complete rewrite from scratch is recommended.**
- **Currently, when submitting jobs via the command line, the entire Rails application is started, and necessary parameters are passed, which is inefficient.** → **Propose creating a dedicated API for job submission.**



