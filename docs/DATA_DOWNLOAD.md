# Climate data acquisition

Level 1 requires no climate downloads. Level 3 uses externally prepared daily
data. The recorded 69-file prepared collection is about 803 GB, separate from
hourly/native archives and workspaces. Installation/tests never download it.

Original source collections:

- CMIP6 climate projections, https://doi.org/10.24381/cds.c866074c
- ERA5 hourly single levels, https://doi.org/10.24381/cds.adbb2d47

Acquire data through the providers' current authorized interfaces with your own
credentials. For each `data/metadata/source_manifest.csv` row select the recorded
model, historical/SSP585 experiment, member, variable, grid and source version.
The manifest's nested source identities retain known version/member identifiers;
missing identifiers remain unknown. Provider access may require registration and
accepting dataset/model terms. Never place API credentials in this repository.

The original request IDs/aggregation receipts are incomplete. Current CDS
downloads are not asserted to reproduce original prepared bytes. BCC historical
r1i1p1f1 v20181126 has an unresolved successor replacement/unpublication chain;
the final used pairing is member-coherent r1i1p1f1. ACCESS historical v20191108
was used, rather than advertised but unavailable v20191118. Do not relabel or
silently replace these versions. Provider availability and permission for every
exact historical version have not been established by this software release.

For **exact-input reproduction**, the corrected prepared files must be obtained
through an author-approved persistent access arrangement, which has not yet been
supplied. The two exact supports are included in `data/supports/`. The upstream collection DOIs alone are not that route.
For **replication from new downloads**, record all actually acquired versions, transformations and
checksums, rebuild supports, and compare resulting statistics without asserting
exact recovery of the article's percentages.

No full climate dataset is embedded in the archive. Public software and compact
results can be archived separately from the external data; that separation must
be described in the final article availability statement with real identifiers.

## What the compact software archive distributes

Included: derived compact CSV/RDS fields and curves, accepted graphics/tables,
histogram supports, source-identity metadata and checksums. Not included: the original prepared
climate collection and large workspaces. The two exact production supports are
included in `data/supports/` (path-redacted copies with unchanged arrays).
Public upstream access routes do not establish exact reconstruction of historical
prepared products: missing acquisition/processing receipts, unavailable original
versions and effective-unit qualifications remain documented. Storage burden,
absence of an exact original copy, missing receipts and source-right restrictions
are separate issues. This archive does not assert a blanket legal prohibition
on sharing prepared data; specific rights must be checked against each provider.
