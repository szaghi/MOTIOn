---
title: stringifor_unicode_tables
---

# stringifor_unicode_tables

> Unicode property tables of StringiFor, generated: do not edit.

**Source**: `src/third_party/StringiFor/src/lib/stringifor_unicode_tables.F90`

**Dependencies**

```mermaid
graph LR
  stringifor_unicode_tables["stringifor_unicode_tables"] --> penf["penf"]
```

## Contents

- [in_ranges](#in-ranges)

## Variables

| Name | Type | Attributes | Description |
|------|------|------------|-------------|
| `UCD_VERSION` | character(len=*) | parameter | Version of the Unicode Character Database. |
| `WIDE_RANGES` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter |  |
| `ZERO_WIDTH_RANGES` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter |  |
| `WHITE_SPACE_RANGES` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | parameter |  |

## Functions

### in_ranges

Return true if the codepoint is in one of the ranges (binary search).

**Attributes**: pure

**Returns**: `logical`

```fortran
function in_ranges(codepoint, ranges) result(found)
```

**Arguments**

| Name | Type | Intent | Attributes | Description |
|------|------|--------|------------|-------------|
| `codepoint` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Codepoint. |
| `ranges` | integer(kind=[I4P](/api/src/third_party/PENF/src/lib/penf_global_parameters_variables)) | in |  | Sorted disjoint ranges, one column [first, last] per range. |

**Call graph**

```mermaid
flowchart TD
  utf8_space_len["utf8_space_len"] --> in_ranges["in_ranges"]
  utf8_unit_width["utf8_unit_width"] --> in_ranges["in_ranges"]
  style in_ranges fill:#3e63dd,stroke:#99b,stroke-width:2px
```
