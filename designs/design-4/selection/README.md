# Candidate-selection evidence

`selected.json.gz`, `selected.tsv` and `selected.fasta` record the 16-region main research selection. `reserve.json.gz` and `reserve.tsv` retain eight additional candidates. Unchanged 25-aa regions shared with Design 2 and Design 3 retain their upstream `region_id`s; shortening `KCTD19` (`event_id` `mn_chr16_67289592_C_T`) from 25 aa to 23 aa assigns a distinct window `region_id` (`REG_6c716df1b65d721ee904` versus Design 2 and Design 3's `REG_9dde223e9bf58444173b`). These are experimental hypotheses, not clinically approved targets.

In `mutation-audit.tsv`, `selected` and `encoded` mean different things and are both retained. `selected` records only whether the event passed the three exploratory gates (paired-normal somatic support, unstranded mutant RNA context, and model nomination), which admits it to the 24-candidate exploratory pool; it is not an encoding decision. `encoded` records whether that event's window is in the panel written to `selected.tsv` for this run. The two differ: `SERPINA1` (`event_id` `mn_chr14_94383173_G_T`) passes the gates but is held in reserve rather than encoded, and `KCTD19` is encoded in the 16-region main panel but not in the 15-region alternative.

RNA counts are read-name fragments under the stated filters. KCTD19 has no coding-strand mutant fragment; FES and SYT12 each have one. Somatic DNA support does not establish expression, presentation or clonality. TJP2 retains transcript/readthrough annotation uncertainty.

Historical fields referring to an SG identifier or an RNA-confirmation watch express an author-selected research hold. The corresponding [research questions](../safety/research-questions.tsv) are not validated assay acceptance criteria, a regulatory protocol or authority to manufacture or administer a product.

The separate 15-region alternative excludes KCTD19 and reoptimizes both the [481-aa protein](../strand-strict-15/construct/protein/cassette.fasta) and its [known-segment PopRNA/GenBank nucleotide encoding](../strand-strict-15/construct/mrna/mrna-segments.fasta). Its clinical benefit and safety are unestablished.
