from incremental_slam.config import RunConfig
from incremental_slam.main import run_selective_demo


def main() -> None:
    cfg = RunConfig(dataset='MIT', batch_size=1, max_gn_iter=10, linear_backend='auto')
    g_final, hist = run_selective_demo(cfg)
    print('Done.')
    print(f'Final nodes: {len(g_final.id_lookup)}')
    print(f'Final edges: {len(g_final.edges)}')
    print(f'Last dx inf-norm: {hist.dx_inf_norm[-1] if hist.dx_inf_norm else 0.0}')
    print(f'Total increments: {len(hist.num_edges)}')

if __name__ == '__main__':
    main()
