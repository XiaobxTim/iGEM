import numpy as np
import pandas as pd

from fitting.fit_module1 import fit_module1_two_stage
from fitting.simulate_module1 import simulate_module1
from fitting.report_module1 import save_module1_report


def main():
    # Example synthetic observation data
    # df = pd.read_csv("blood_data.csv")

    # observations = {
    #     "A_blood": {
    #         "t": df["time_h"].to_numpy(dtype=float),
    #         "y": df["A_blood"].to_numpy(dtype=float),
    #         "weight": 1.0,
    #     }
    # }

    observations = {
        "A_blood": {
            "t": np.array([1, 4, 8, 24, 48, 72], dtype=float),
            "y": np.array([0.03, 0.10, 0.16, 0.24, 0.21, 0.14], dtype=float),
            "weight": 1.0,
        }
    }

    dose = 1.0

    fit_result = fit_module1_two_stage(
        dose=dose,
        observations=observations,
        loss_type="log_sse",
        use_soft_penalty=False,
    )

    # fit_result = fit_module1_two_stage(
    #     dose=dose,
    #     observations=observations,
    #     loss_type="log_sse",
    #     use_soft_penalty=True,
    #     penalty_weight=10.0,
    # )

    print("Best parameters:")
    for k, v in fit_result["best_params"].items():
        print(f"  {k}: {v:.6f}")
    print(f"Best loss: {fit_result['best_loss']:.6f}")

    t_eval = np.linspace(0, 72, 400)
    sim = simulate_module1(
        params=fit_result["best_params"],
        dose=dose,
        t_eval=t_eval,
    )

    idx = sim["idx"]
    print("Final blood value:", sim["y"][idx["A_blood"], -1])

    saved = save_module1_report(
        fit_result=fit_result,
        dose=dose,
        observations=observations,
        output_dir="outputs/module1_fit_report",
        t_end=72.0,
    )

    print("Saved outputs:")
    for name, path in saved.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()