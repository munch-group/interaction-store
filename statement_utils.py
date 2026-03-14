def load_store(path=None):
    """Load statements from JSON file, return list of INDRA Statement objects."""
    if not path.exists():
        return []
    with open(path) as f:
        return stmts_from_json(json.load(f))


def save_store(stmts, path=None):
    """Serialise statement list to JSON file."""
    with open(path, 'w') as f:
        json.dump(stmts_to_json(stmts), f, indent=2)
    print(f'Saved {len(stmts)} statements → {path}')


def add_statements(new_stmts, path=None):
    """Append new statements to the persistent store."""
    existing = load_store(path)
    combined = existing + new_stmts
    save_store(combined, path)
    return combined


def ev(text, pmid=None, doi=None, context=None, hypothesis=False, direct=True):
    """Convenience: create a manual Evidence object.

    Args:
        text:       Mechanistic reasoning or supporting sentence.
        pmid:       PubMed ID string, e.g. '23640493'.
        doi:        DOI string, e.g. '10.1161/ATVBAHA.113.301522'.
        context:    Module name or thematic tag.
        hypothesis: True if this is speculative / not yet confirmed.
        direct:     False if interaction is pathway-level, not protein-protein.
    """
    text_refs = {}
    if pmid: text_refs['PMID'] = str(pmid)
    if doi:  text_refs['DOI']  = doi
    return Evidence(
        text=text,
        source_api='manual',
        pmid=pmid,
        text_refs=text_refs or None,
        annotations={
            'date':       str(date.today()),
            'context':    context or 'exploratory',
            'directness': 'direct' if direct else 'indirect',
        },
        epistemics={'hypothesis': hypothesis},
    )


def ag(name, hgnc=None, up=None):
    """Convenience: create an Agent with optional db_refs."""
    db_refs = {}
    if hgnc:
        db_refs['HGNC'] = str(hgnc)
    if up:
        db_refs['UP'] = up
    return Agent(name, db_refs=db_refs or None)


def summarise(path=None):
    """Print a readable summary of the store."""
    stmts = load_store(path)
    print(f'{len(stmts)} statements in store\n')
    by_type = {}
    for s in stmts:
        t = type(s).__name__
        by_type.setdefault(t, []).append(s)
    for t, ss in sorted(by_type.items()):
        print(f'  {t}: {len(ss)}')
    print()
    for s in stmts:
        subj = s.agent_list()[0].name if s.agent_list() else '?'
        obj  = s.agent_list()[-1].name if len(s.agent_list()) > 1 else ''
        txt  = s.evidence[0].text[:70] if s.evidence else ''
        print(f'  [{type(s).__name__:20s}]  {subj:12s} → {obj:12s}  "{txt}"')
