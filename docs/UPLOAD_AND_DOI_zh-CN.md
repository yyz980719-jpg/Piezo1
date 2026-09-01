# GitHub + Zenodo 上传与 DOI 注册

## 推荐路线：GitHub release 自动归档到 Zenodo

1. 在 GitHub 新建公开仓库，不要自动添加 README、LICENSE 或 `.gitignore`。
2. 将本目录的内容上传到仓库根目录；不要只上传外层压缩包，也不要形成双层同名目录。
3. 推送 `main` 分支，但暂时不要创建 release。
4. 登录 Zenodo，连接 GitHub；在 Zenodo 的 GitHub 页面点击 **Sync now**，找到该仓库并启用。
5. 回到 GitHub，以已经存在的标签 `v2.0.0` 创建 release；标题建议为
   `PIEZO1 OA reproducible research compendium v2.0.0`。
6. 等待 Zenodo 自动摄取，核对作者、ORCID、单位、版本、MIT 许可和文件清单后发布记录。
7. 在无登录状态下测试 DOI 和文件下载。
8. 将 DOI 补入 `CITATION.cff`、GitHub README 和论文 Code Availability；不要移动或重建已归档的 `v2.0.0` 标签。

官方说明：

- Zenodo 启用 GitHub 仓库：https://help.zenodo.org/docs/github/enable-repository/
- GitHub `CITATION.cff`：https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files

## 备选路线：Zenodo 手动上传

上传根目录提供的发布 ZIP。在 Zenodo 草稿中先核对元数据并保留 DOI，再发布。
若先保留 DOI，可在最终上传前把 DOI 写入 CFF；任何重新生成的 ZIP 都必须重新计算 SHA-256。

## 发布前必须核对

- Git 标签与拟发布版本均为 `v2.0.0`。
- Git 工作区干净，标签指向最终提交。
- 压缩包不含 `.git/`、原始第三方大数据、账号令牌或本地绝对路径。
- `CITATION.cff` 和 `.zenodo.json` 中作者、ORCID、单位一致。
- Zenodo 文件列表中能看到 README、LICENSE、CFF、代码、`renv.lock`、reference 和 manuscript。
- DOI 真正解析后，再在论文中声称代码公开可用。
