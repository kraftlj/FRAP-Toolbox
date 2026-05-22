# Data Availability And Fixture Manifest

The source repository is hosted at
[https://github.com/kraftlj/FRAP-Toolbox](https://github.com/kraftlj/FRAP-Toolbox).
The old project website is no longer used as an active download location in the
documentation or software citation.

GitHub is not the right home for the full raw microscopy fixture archive: the
current local `test-data/` folder is about 1.2 GB, several files are close to
100 MB, and GitHub warns on regular Git files above 50 MiB while blocking files
above 100 MiB. The repository therefore tracks one raw sample for smoke tests
and keeps full parity fixtures external.

## Tracked Sample

| Path | Size bytes | Size MiB | SHA-256 | Use |
| --- | ---: | ---: | --- | --- |
| `sample-data/Diffusion/Venus_Cytoplasm_1.lsm` | 31325886 | 29.87 | `4f23892b26d19e9611b5ec76e8675cc66c8c7962a043d84372ba535545f1f804` | Canonical diffusion raw-image smoke sample and first original guide example. |

## Zenodo Record

The full fixture archive is published as a Zenodo dataset:

- Title: `FRAP-Toolbox Legacy User Guide Test Data`
- Creators: Lewis J. Kraft; Jacob Dowler; Charles A. Day; Minchul Kang; Anne K. Kenworthy
- Publication date: 2014-06-14
- License: `CC-BY-4.0`
- Version DOI: [https://doi.org/10.5281/zenodo.20344310](https://doi.org/10.5281/zenodo.20344310)
- All-versions DOI: [https://doi.org/10.5281/zenodo.20344309](https://doi.org/10.5281/zenodo.20344309)
- Record page: [https://zenodo.org/records/20344310](https://zenodo.org/records/20344310)
- Archive name: `frap-toolbox-legacy-user-guide-test-data.zip`
- Archive SHA-256: `3c17d0453bd325d3032859af75e90146008a7612a7bcc6fa9b92d90a9d046ddc`
- Archive size: 845.4 MB on Zenodo; about 1.34 GB uncompressed

## Restoring The Full Archive

Download `frap-toolbox-legacy-user-guide-test-data.zip` from the Zenodo record,
then unpack it at the repository root so paths match the legacy guide and the
optional parity tests:

```text
test-data/
  Userguide.pdf
  Diffusion/
  Reaction 1/
  Reaction 2/
```

The `test-data/` directory remains ignored by git. Tests that require it must
skip when the directory is absent.

## Full Archive Manifest

The table below records the expected local `test-data/` contents from the
current legacy archive. `.DS_Store` and other local metadata files are excluded.

| Source path | Size bytes | Size MiB | SHA-256 | Use |
| --- | ---: | ---: | --- | --- |
| `test-data/Diffusion/Venus-ATG5_Cytoplasm_Diffusion_FRAP_datasets.txt` | 622086 | 0.59 | `20c0f23d57740b4829d44c5e38483932a25f04676fd2620b405111a719d2456d` | Diffusion exported FRAP vectors for model/parity tests. |
| `test-data/Diffusion/Venus-ATG5_Cytoplasm_Diffusion_Fit_Parameters.txt` | 1416 | 0.00 | `2fb6b0f5376e31c8c8a96d495bf14f070f31d97acaa07f608cb81b4d34bce6b8` | Diffusion guide parameter-table contract. |
| `test-data/Diffusion/Venus-ATG5_Cytoplasm_Diffusion_Postbleach_profiles.txt` | 3359216 | 3.20 | `e0bf003eacb9ba9ea05a7352cd98b0a4c36dbaa1c73b192c44dd3c076efd822c` | Diffusion exported post-bleach profiles for model/parity tests. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_1.lsm` | 31325896 | 29.87 | `1594c14f2a1c2e90139c22f476ff2e80528bce9d74cf653c9c4307760ee951c6` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_10.lsm` | 31325898 | 29.87 | `8076e48271d40c74833b8abd4af09fdba95d79a3355e101f03e62353de6b012d` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_11.lsm` | 31325898 | 29.87 | `e70d18b08e30128944e0292f2438e705e72c9db23e0850dc372f328467cdd4ae` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_12.lsm` | 31325898 | 29.87 | `4278bf792531178a562c1fdbaeb6940c121898cab153f2b25670c514772020ba` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_13.lsm` | 31325898 | 29.87 | `00264b4b9f804a2d5d175abac306b9628e9d6424673f258eeec7faa52c69996a` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_14.lsm` | 31325898 | 29.87 | `fc1c588f34a9605c071245b344ff1dc86c69b17cdc2a2102af3b3865bb1f82f2` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_15.lsm` | 31325898 | 29.87 | `a50b7c7d63baadc2840e59cd61841515d497ef049c4366da2839e2e3b6bfd70b` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_2.lsm` | 31325896 | 29.87 | `fb9d69fb4c9e7dc27e2cd6a7b8dcfb3c58a5808b45cf29847e4d027950c2437a` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_3.lsm` | 31325896 | 29.87 | `a89f5b36982396a7cd9aee7701b1534a60177043b629a506cbb699a10008305f` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_4.lsm` | 31325896 | 29.87 | `3a06364b5f0f9439602f7c532c2615b321adf97493aae3de16667c76841ff26b` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_5.lsm` | 31325896 | 29.87 | `5d45078923f13363fd35762fd12780f3fb42ef6d56be6fa4fa93ea6d7d4f75b2` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_6.lsm` | 31325896 | 29.87 | `1f2341068463951199eb08256c2dddecf75d823b28eefdb6382091774028a506` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_7.lsm` | 31325896 | 29.87 | `c0df470ec09e950c6bb4dbfadf2bf370de8fdb4bdd31bab2be4308504c22e822` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_8.lsm` | 31325896 | 29.87 | `7b84d4c4e6899e4cb071ffceffcf77237787261943700f840d459174c14f2e48` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus-Atg5_Cytoplasm_9.lsm` | 31325896 | 29.87 | `94de62ab59b65de69fa6f251afb5102ab7eb3712d244e5822f187a7221396242` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_1.lsm` | 31325886 | 29.87 | `4f23892b26d19e9611b5ec76e8675cc66c8c7962a043d84372ba535545f1f804` | Canonical diffusion raw-image smoke/parity sample; copied to sample-data for tracked demos. |
| `test-data/Diffusion/Venus_Cytoplasm_10.lsm` | 31325888 | 29.87 | `aef3c8e5661853315005fe08c8a4af127c5222af37bc218e419fed86a0a8c148` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_11.lsm` | 31325888 | 29.87 | `d03bb3866633897b02253c0619198aa84bc14f4ce00958204945135f0e500f45` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_12.lsm` | 30405364 | 29.00 | `787f2a38b2a12f3cd220bedd46d71a98ae3d2a49672d6f240f3c541593c5e459` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_13.lsm` | 31325888 | 29.87 | `16e2a10e0b145631fd9b6612619ef9c97149453bb0088addd106a13c45b9ece2` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_14.lsm` | 31325888 | 29.87 | `0471b618cbb59481a9de1bc6dfb12d8f919941fcdfd83e8100d237e275f7c086` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_15.lsm` | 31325888 | 29.87 | `a5fb3123eaec81f7dac3b484226d3917bc235777906e4b28c5dfcc475ca4cf9f` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_2.lsm` | 31325886 | 29.87 | `39fa3ebad1963885562fe1df7637dec7d034ed9297a6b6dfd189f4c40557707f` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_3.lsm` | 31325886 | 29.87 | `e9f9712d567d3a52275403bcdfcf088c048328860121f8c7593ebe70e32bbd10` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_4.lsm` | 31325886 | 29.87 | `08518bb7e8f2bdc39ced5675875200e48dd9e14c4fa199a0bc8d350c450550a9` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_5.lsm` | 31325886 | 29.87 | `0143ab3968f5c6e43a4e3f620b5b6b6dd817363c903af09d2ab4a66c8ff183c3` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_6.lsm` | 31325886 | 29.87 | `710056f2b3db7a81cbb4688e8da8e1ca00f7b588a91a6ba3e6ed0160c22f07dd` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_7.lsm` | 31325886 | 29.87 | `79aef8834d42eb8fc80a91e0277f086be02abc00cfdcd10371a6de0982b127cf` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_8.lsm` | 31325886 | 29.87 | `b3e18337de9e0b8311732528445bdbb04a87e272bc3e6270adaf9a852a32e019` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_9.lsm` | 31325886 | 29.87 | `663adf1b63d164a79b6b8e897ea1010a35065923feb8f5224d8930bdd95e483b` | Full diffusion raw-image parity fixture. |
| `test-data/Diffusion/Venus_Cytoplasm_Diffusion_FRAP_datasets.txt` | 626072 | 0.60 | `45a0c6f475690896e74b72edd034338ade01d22e347929be3ffded52c998a80f` | Diffusion exported FRAP vectors for model/parity tests. |
| `test-data/Diffusion/Venus_Cytoplasm_Diffusion_Fit_Parameters.txt` | 1340 | 0.00 | `63a6bc6d4a9c540219398be8ce02886c5b58c63e01f9294be6c8ae3b0b8f296d` | Diffusion guide parameter-table contract. |
| `test-data/Diffusion/Venus_Cytoplasm_Diffusion_Postbleach_profiles.txt` | 3354580 | 3.20 | `6fe7ce9a6eb5d99b43cb703fb57c1d824cfb3fad9490f729cf5f3863c7b069f3` | Diffusion exported post-bleach profiles for model/parity tests. |
| `test-data/Reaction 1/Venus_1001.nd2` | 97882112 | 93.35 | `d189df1d91d3bcc5e44d15879dda5c5bfa412389ced5beef158f555df758a835` | Reaction 1 raw ND2 fixture for loader diagnostics and future ROI parity. |
| `test-data/Reaction 1/Venus_1002.nd2` | 97964032 | 93.43 | `c889e53e886b68ed7b245774ab2f24f43a41c07b9182b3ce8d6a7e1b721da46b` | Reaction 1 raw ND2 fixture for loader diagnostics and future ROI parity. |
| `test-data/Reaction 1/Venus_NCTransport_Reaction_FRAP_datasets.txt` | 34702 | 0.03 | `6e5fb1bb48386f99829dd44f946b1d306f781f81a98716126324830f7df1c8b9` | Reaction 1 exported FRAP vectors for refit parity. |
| `test-data/Reaction 1/Venus_NCTransport_Reaction_Fit_Parameters.txt` | 194 | 0.00 | `bc17f834870c8a48fb28ab6dfa8ce85777c156551368171bf9d2df6a0d2cf1b8` | Reaction 1 guide parameter-table contract. |
| `test-data/Reaction 2/Venus-Atg5_1002.nd2` | 97968128 | 93.43 | `7e6df9bc9f528b70f48ba9dda470add65109fb2c799ac69c958649ed1883cd7e` | Reaction 2 raw ND2 fixture for loader diagnostics and future ROI parity. |
| `test-data/Reaction 2/Venus-Atg5_1003.nd2` | 97964032 | 93.43 | `600092d98af6a097f4a91d9e32af0e2320868725c2292d94609011fb8eb357a9` | Reaction 2 raw ND2 fixture for loader diagnostics and future ROI parity. |
| `test-data/Reaction 2/Venus-Atg5_NCTransport_Reaction2_Fit_Parameters.txt` | 166 | 0.00 | `0bedb9e37bc90f1a70ee701e1e3683723e52b197933054ee796638b03d29668e` | Reaction 2 guide parameter-table contract. |
| `test-data/Userguide.pdf` | 440120 | 0.42 | `9d4b2b02ad8adc6bbaeeffe6b0864629cb189cbf280a52741e8ef439914624fd` | Original PDF source for Markdown guide recreation; parity reference only. |

## Hosting Notes

Relevant limits and policies:

- GitHub regular Git files: warning above 50 MiB, blocked above 100 MiB.
- GitHub recommended repository size: ideally below 1 GB.
- Zenodo default record limit: 50 GB and 100 files, with quota increases by
  request.
- OSF considered alternative: public OSF Storage supports larger projects than
  private storage, but Zenodo gives the cleanest DOI-first archive path for this
  fixture set.
