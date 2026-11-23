"""Training-only network selection and shared progress descriptions."""
import math


def select(report, objective=None):
    """Describe retained weights; evaluation metrics never select them.

    Historical reports contain objective history but may lack an explicit
    selection record. An optional checkpoint objective identifies the exact
    retained minimum, including when a later iteration was interrupted.
    """
    saved=report.get('model_selection')
    if saved is not None and objective is None:
        return dict(saved)
    candidates=[]
    for row in report.get('history',[]):
        value=row.get('objective')
        if value is None or not math.isfinite(value):continue
        if objective is not None and value!=objective:continue
        candidates.append(row)
    if not candidates:return None
    row=min(candidates,key=lambda x:x['objective'])
    return dict(iteration=int(row['iteration']),
                objective=float(row['objective']),
                criterion='minimum mean training-basin objective',
                basin_group='training',period='training',
                loss_id=report.get('config',{}).get('loss'),
                evaluation_used=False)


def record(report,item):
    """Record the network evaluated at this iteration, before its update."""
    selection=select(dict(config=report['config'],history=[item]))
    if selection is None:raise ValueError('Cannot retain a nonfinite objective.')
    report['model_selection']=selection
    return selection


def completed(report):
    """Count completed evaluations rather than an interrupted proposal."""
    return max((int(x['iteration']) for x in report.get('history',[])
                if x.get('iteration_complete') or x.get('iteration_text')
                or x.get('diagnostics')),default=0)


def progress_text(done,total,stage=14,finished=False,status='completed',current=None):
    message=f'({stage:02d}) train shared network ...'
    if finished:
        return message+f' {status} {done:,} iterations.'
    counter=f'iteration {current:,}' if current is not None else f'{done:,}'
    return message+f' [{counter}/{total:,}, {100*done/max(1,total):.1f}% done]'


def selection_text(selection,last=None):
    if not selection:return '      No completed training network was selected.'
    iteration=selection['iteration']
    last=last if last is not None else selection.get('completed_iterations')
    count=f' of {last:,}' if last is not None else ''
    return (f'      Selected iteration {iteration:,}{count}:\n'
            f'    minimum mean training loss = {selection["objective"]:.9g} '
            '(train basins | train period)')
