from __future__ import annotations

import math
from typing import List, Tuple

import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from atr_sim.core.geometry import Polygon, Ray
from atr_sim.core import physics
from atr_sim.utils.constants import REFRACTIVE_INDICES


class Application(tk.Tk):
    def __init__(self) -> None:
        super().__init__()
        self.title("ATR-MIR Ray Path Simulator")
        self.geometry("1100x700")

        self._build_ui()

    # UI construction
    def _build_ui(self) -> None:
        self.columnconfigure(0, weight=0)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)

        # Control panel
        control_frame = ttk.Frame(self, padding=10)
        control_frame.grid(row=0, column=0, sticky="nsw")

        # Plot panel
        plot_frame = ttk.Frame(self, padding=10)
        plot_frame.grid(row=0, column=1, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1)
        plot_frame.columnconfigure(0, weight=1)

        # Shape selection and dimensions
        shape_label = ttk.Label(control_frame, text="Shape:")
        shape_label.grid(row=0, column=0, sticky="w")
        self.shape_var = tk.StringVar(value="Rectangle")
        shape_combo = ttk.Combobox(control_frame, textvariable=self.shape_var, values=["Rectangle", "Trapezoid"], state="readonly", width=14)
        shape_combo.grid(row=0, column=1, sticky="w", pady=2)
        shape_combo.bind("<<ComboboxSelected>>", lambda e: self._on_shape_change())

        # Rectangle dimensions
        self.rect_width_var = tk.StringVar(value="50.0")
        self.rect_height_var = tk.StringVar(value="20.0")
        ttk.Label(control_frame, text="Rect Width:").grid(row=1, column=0, sticky="w")
        self.rect_width_entry = ttk.Entry(control_frame, textvariable=self.rect_width_var, width=16)
        self.rect_width_entry.grid(row=1, column=1, sticky="w", pady=2)
        ttk.Label(control_frame, text="Rect Height:").grid(row=2, column=0, sticky="w")
        self.rect_height_entry = ttk.Entry(control_frame, textvariable=self.rect_height_var, width=16)
        self.rect_height_entry.grid(row=2, column=1, sticky="w", pady=2)

        # Trapezoid dimensions (with constraint)
        self.trap_base_var = tk.StringVar(value="30.0")
        self.trap_top_var = tk.StringVar(value="10.0")
        self.trap_constrain_var = tk.StringVar(value="Height")  # Height or Angle
        self.trap_height_var = tk.StringVar(value="20.0")
        self.trap_angle_deg_var = tk.StringVar(value="45.0")
        ttk.Label(control_frame, text="Trap Base L:").grid(row=3, column=0, sticky="w")
        self.trap_base_entry = ttk.Entry(control_frame, textvariable=self.trap_base_var, width=16)
        self.trap_base_entry.grid(row=3, column=1, sticky="w", pady=2)
        ttk.Label(control_frame, text="Trap Top L:").grid(row=4, column=0, sticky="w")
        self.trap_top_entry = ttk.Entry(control_frame, textvariable=self.trap_top_var, width=16)
        self.trap_top_entry.grid(row=4, column=1, sticky="w", pady=2)
        ttk.Label(control_frame, text="Constrain by:").grid(row=5, column=0, sticky="w")
        self.trap_constrain_combo = ttk.Combobox(control_frame, textvariable=self.trap_constrain_var, values=["Height", "Angle"], state="readonly", width=14)
        self.trap_constrain_combo.grid(row=5, column=1, sticky="w", pady=2)
        self.trap_constrain_combo.bind("<<ComboboxSelected>>", lambda e: self._on_trap_constraint_change())
        ttk.Label(control_frame, text="Trap Height:").grid(row=6, column=0, sticky="w")
        self.trap_height_entry = ttk.Entry(control_frame, textvariable=self.trap_height_var, width=16)
        self.trap_height_entry.grid(row=6, column=1, sticky="w", pady=2)
        ttk.Label(control_frame, text="Top Angle (deg):").grid(row=7, column=0, sticky="w")
        self.trap_angle_entry = ttk.Entry(control_frame, textvariable=self.trap_angle_deg_var, width=16)
        self.trap_angle_entry.grid(row=7, column=1, sticky="w", pady=2)

        ttk.Separator(control_frame, orient="horizontal").grid(row=8, column=0, columnspan=2, sticky="ew", pady=6)

        # Material selection
        ttk.Label(control_frame, text="Crystal Material:").grid(row=9, column=0, sticky="w")
        self.material_var = tk.StringVar(value="ZnSe")
        material_options = list(REFRACTIVE_INDICES.keys())
        material_combo = ttk.Combobox(control_frame, textvariable=self.material_var, values=material_options, state="readonly", width=14)
        material_combo.grid(row=9, column=1, sticky="w", pady=2)
        material_combo.bind("<<ComboboxSelected>>", lambda e: self._on_material_change())

        ttk.Label(control_frame, text="Crystal n:").grid(row=10, column=0, sticky="w")
        self.material_n_var = tk.StringVar(value=str(REFRACTIVE_INDICES[self.material_var.get()]))
        self.material_n_entry = ttk.Entry(control_frame, textvariable=self.material_n_var, width=16)
        self.material_n_entry.grid(row=10, column=1, sticky="w", pady=2)

        ttk.Separator(control_frame, orient="horizontal").grid(row=11, column=0, columnspan=2, sticky="ew", pady=6)

        # Boundary media indices (4 sides)
        self.boundary_vars: List[tk.StringVar] = [tk.StringVar(value="1.0") for _ in range(4)]
        for i in range(4):
            ttk.Label(control_frame, text=f"Boundary {i} n:").grid(row=12 + i, column=0, sticky="w")
            ttk.Entry(control_frame, textvariable=self.boundary_vars[i], width=16).grid(row=12 + i, column=1, sticky="w", pady=2)

        ttk.Separator(control_frame, orient="horizontal").grid(row=16, column=0, columnspan=2, sticky="ew", pady=6)

        # Ray parameters
        ttk.Label(control_frame, text="Ray Origin X:").grid(row=17, column=0, sticky="w")
        self.ray_x_var = tk.StringVar(value="5.0")
        ttk.Entry(control_frame, textvariable=self.ray_x_var, width=16).grid(row=17, column=1, sticky="w", pady=2)

        ttk.Label(control_frame, text="Ray Origin Y:").grid(row=18, column=0, sticky="w")
        self.ray_y_var = tk.StringVar(value="-5.0")
        ttk.Entry(control_frame, textvariable=self.ray_y_var, width=16).grid(row=18, column=1, sticky="w", pady=2)

        ttk.Label(control_frame, text="Ray Angle (deg):").grid(row=19, column=0, sticky="w")
        self.ray_angle_deg_var = tk.StringVar(value="90.0")
        ttk.Entry(control_frame, textvariable=self.ray_angle_deg_var, width=16).grid(row=19, column=1, sticky="w", pady=2)

        ttk.Separator(control_frame, orient="horizontal").grid(row=20, column=0, columnspan=2, sticky="ew", pady=6)

        # Actions
        run_btn = ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation)
        run_btn.grid(row=21, column=0, columnspan=2, sticky="ew", pady=4)

        clear_btn = ttk.Button(control_frame, text="Clear", command=self.clear_plot)
        clear_btn.grid(row=22, column=0, columnspan=2, sticky="ew", pady=2)

        # Toggle interactions visibility (to reduce side panel width when hidden)
        self.show_info_var = tk.BooleanVar(value=True)
        self.toggle_info_btn = ttk.Checkbutton(
            control_frame,
            text="Show interactions",
            variable=self.show_info_var,
            command=self._toggle_info,
        )
        self.toggle_info_btn.grid(row=23, column=0, columnspan=2, sticky="w", pady=(2, 2))

        # Interaction info section
        self.info_frame = ttk.LabelFrame(control_frame, text="Interactions", padding=(6, 6))
        self.info_frame.grid(row=24, column=0, columnspan=2, sticky="nsew", pady=6)
        self.info_frame.columnconfigure(0, weight=1)

        # Compact treeview styling
        style = ttk.Style(self)
        try:
            style.configure("Compact.Treeview", rowheight=16, font=("Segoe UI", 8))
            style.configure("Compact.Treeview.Heading", font=("Segoe UI", 8, "bold"))
        except Exception:
            pass

        columns = ("hit", "boundary", "mode", "x", "y", "theta_i", "critical", "n_pair", "dir")
        self.info_tree = ttk.Treeview(self.info_frame, columns=columns, show="headings", height=5, style="Compact.Treeview")
        # Short headings to reduce width
        self.info_tree.heading("hit", text="#")
        self.info_tree.heading("boundary", text="B")
        self.info_tree.heading("mode", text="Mode")
        self.info_tree.heading("x", text="X")
        self.info_tree.heading("y", text="Y")
        self.info_tree.heading("theta_i", text="θi")
        self.info_tree.heading("critical", text="θc")
        self.info_tree.heading("n_pair", text="n")
        self.info_tree.heading("dir", text="Dir")
        # Narrower column widths
        self.info_tree.column("hit", width=28, anchor="center")
        self.info_tree.column("boundary", width=70, anchor="w")
        self.info_tree.column("mode", width=46, anchor="center")
        self.info_tree.column("x", width=65, anchor="e")
        self.info_tree.column("y", width=65, anchor="e")
        self.info_tree.column("theta_i", width=70, anchor="e")
        self.info_tree.column("critical", width=70, anchor="e")
        self.info_tree.column("n_pair", width=90, anchor="center")
        self.info_tree.column("dir", width=80, anchor="center")
        self.info_tree.grid(row=0, column=0, sticky="nsew")

        # Optional scrollbars without expanding width
        yscroll = ttk.Scrollbar(self.info_frame, orient="vertical", command=self.info_tree.yview)
        self.info_tree.configure(yscrollcommand=yscroll.set)
        yscroll.grid(row=0, column=1, sticky="ns")

        # Boundary mapping helper label (shortened)
        self.boundary_map_label = ttk.Label(self.info_frame, text="", justify="left")
        self.boundary_map_label.grid(row=1, column=0, sticky="w", pady=(6, 0))

        # Plot setup
        self.figure = Figure(figsize=(6, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        # Initial view; will be auto-zoomed after plotting
        self.ax.set_xlim(-50, 150)
        self.ax.set_ylim(-50, 150)
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")
        self.toolbar = NavigationToolbar2Tk(self.canvas, plot_frame, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.grid(row=1, column=0, sticky="ew")

        self._on_shape_change()  # set initial visibility

    def _on_shape_change(self) -> None:
        shape = self.shape_var.get()
        # Toggle visibility of rectangle vs trapezoid fields
        rect_widgets = [self.rect_width_entry, self.rect_height_entry]
        trap_widgets = [
            self.trap_base_entry,
            self.trap_top_entry,
            self.trap_constrain_combo,
            self.trap_height_entry,
            self.trap_angle_entry,
        ]

        if shape == "Rectangle":
            for w in rect_widgets:
                w.configure(state="normal")
            for w in trap_widgets:
                w.configure(state="disabled")
        else:
            for w in rect_widgets:
                w.configure(state="disabled")
            for w in trap_widgets:
                w.configure(state="normal")
            # Apply constraint toggle for trapezoid
            self._on_trap_constraint_change()

    def _on_trap_constraint_change(self) -> None:
        # Only one of height or angle should be editable at a time
        mode = self.trap_constrain_var.get()
        if mode == "Height":
            self.trap_height_entry.configure(state="normal")
            self.trap_angle_entry.configure(state="disabled")
        else:
            self.trap_height_entry.configure(state="disabled")
            self.trap_angle_entry.configure(state="normal")

    def _on_material_change(self) -> None:
        key = self.material_var.get()
        if key in REFRACTIVE_INDICES:
            self.material_n_var.set(str(REFRACTIVE_INDICES[key]))

    def clear_plot(self) -> None:
        self.ax.clear()
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        # Neutral default view before next simulation
        self.ax.set_xlim(-50, 150)
        self.ax.set_ylim(-50, 150)
        self.canvas.draw()
        # Clear interaction info
        if hasattr(self, "info_tree"):
            for item in self.info_tree.get_children():
                self.info_tree.delete(item)
        if hasattr(self, "boundary_map_label"):
            self.boundary_map_label.configure(text="")

    def _compute_vertices(self) -> List[Tuple[float, float]]:
        shape = self.shape_var.get()
        if shape == "Rectangle":
            width = float(self.rect_width_var.get())
            height = float(self.rect_height_var.get())
            if width <= 0 or height <= 0:
                raise ValueError("Rectangle dimensions must be positive")
            return [(0.0, 0.0), (width, 0.0), (width, height), (0.0, height)]
        else:
            base_L = float(self.trap_base_var.get())
            top_L = float(self.trap_top_var.get())
            if base_L <= 0:
                raise ValueError("Trapezoid base must be positive")
            if not (0.0 < top_L < base_L):
                raise ValueError("Trap Top L must be between 0 and base length")

            mode = self.trap_constrain_var.get()
            if mode == "Height":
                h = float(self.trap_height_var.get())
                if h <= 0:
                    raise ValueError("Trapezoid height must be positive")
                offset = 0.5 * (base_L - top_L)
            else:
                angle_deg = float(self.trap_angle_deg_var.get())
                angle_rad = math.radians(angle_deg)
                t = math.tan(angle_rad)
                if t <= 0:
                    raise ValueError("Angle must be > 0 deg")
                offset = 0.5 * (base_L - top_L)
                h = offset / t
                if h <= 0:
                    raise ValueError("Computed height is non-positive; check inputs")

            # CCW vertices
            return [
                (0.0, 0.0),
                (base_L, 0.0),
                (base_L - offset, h),
                (offset, h),
            ]

    def get_params_from_gui(self):
        # Geometry
        vertices = self._compute_vertices()
        material_n = float(self.material_n_var.get())
        boundary_ns = [float(v.get()) for v in self.boundary_vars]

        # Ray
        rx = float(self.ray_x_var.get())
        ry = float(self.ray_y_var.get())
        theta_deg = float(self.ray_angle_deg_var.get())
        theta = math.radians(theta_deg)
        ray_origin = np.array([rx, ry], dtype=float)
        ray_dir = np.array([math.cos(theta), math.sin(theta)], dtype=float)

        return vertices, material_n, boundary_ns, ray_origin, ray_dir

    def run_simulation(self) -> None:
        try:
            vertices, material_n, boundary_ns, ray_origin, ray_dir = self.get_params_from_gui()
        except ValueError as e:
            messagebox.showerror("Invalid input", str(e))
            return

        polygon = Polygon(vertices, material_n, default_boundary_n=boundary_ns[0])
        # Ensure we set per-boundary values using UI→polygon index mapping
        index_map = self._boundary_index_map(polygon)
        for ui_idx, poly_idx in enumerate(index_map):
            polygon.set_boundary_medium(poly_idx, boundary_ns[ui_idx % len(boundary_ns)])

        # Prepare interaction info table
        for item in getattr(self, "info_tree", []).get_children() if hasattr(self, "info_tree") else []:
            self.info_tree.delete(item)
        self._update_boundary_mapping(polygon)

        current_ray = Ray(ray_origin, ray_dir)
        ray_path: List[np.ndarray] = [current_ray.origin.copy()]
        MAX_BOUNCES = 50
        EXIT_DISTANCE = 10.0
        EPS = 1e-7

        hit_index = 1
        for _ in range(MAX_BOUNCES):
            hit = physics.find_closest_intersection(current_ray, polygon, epsilon=1e-9)
            if hit is None:
                # No more intersections found; extend ray forward a bit and stop
                final_point = current_ray.origin + current_ray.direction * EXIT_DISTANCE
                ray_path.append(final_point)
                break

            new_start = hit["intersection_point"]
            ray_path.append(new_start)
            normal = hit["normal"]
            segment_idx = hit["segment_index"]

            n1 = polygon.material_n
            n2 = polygon.boundary_media_n[segment_idx]

            # Debug incident/critical angles
            dbg = physics.compute_incident_debug(current_ray.direction, normal, n1, n2)
            angle_i = dbg["angle_i_deg"]
            crit = dbg["critical_deg"]
            if dbg["inside_to_outside"] == 1.0:
                direction_txt = "inside→out"
            else:
                direction_txt = "outside→in"

            new_direction, mode, inside_to_outside = physics.calculate_interaction(current_ray.direction, normal, n1, n2)

            # Override direction text for TIR cases
            if mode == "reflect":
                direction_txt = "inside -> inside"

            # Append info row
            if hasattr(self, "info_tree"):
                boundary_name = self._boundary_name(segment_idx)
                crit_txt = "—" if crit < 0 else f"{crit:.2f}"
                self.info_tree.insert("", "end", values=(
                    hit_index,
                    f"{segment_idx} {boundary_name}",
                    "TIR" if mode == "reflect" else mode,
                    f"{new_start[0]:.3f}",
                    f"{new_start[1]:.3f}",
                    f"{angle_i:.2f}",
                    crit_txt,
                    f"{dbg['n_incident']:.3f}→{dbg['n_transmitted']:.3f}",
                    direction_txt,
                ))
                hit_index += 1

            if mode == "refract" and inside_to_outside:
                # Ray exits the crystal; show a short exit segment and stop
                exit_point = new_start + new_direction * EXIT_DISTANCE
                ray_path.append(exit_point)
                break

            # Continue tracing (reflection or internal refraction)
            # Offset slightly along new direction to avoid self-intersection with the same edge
            current_ray = Ray(new_start + new_direction * EPS, new_direction)

        self.plot_results(polygon, ray_path)

    def plot_results(self, polygon: Polygon, ray_path: List[np.ndarray]) -> None:
        self.ax.clear()
        self.ax.set_aspect('equal', adjustable='box')
        # Plot polygon
        xs = [p[0] for p in polygon.vertices] + [polygon.vertices[0][0]]
        ys = [p[1] for p in polygon.vertices] + [polygon.vertices[0][1]]
        self.ax.plot(xs, ys, 'k-', linewidth=2, label='Crystal')

        # Plot ray path
        if len(ray_path) >= 2:
            rxs = [p[0] for p in ray_path]
            rys = [p[1] for p in ray_path]
            self.ax.plot(rxs, rys, 'r-', linewidth=1.5, label='Ray path')
            self.ax.scatter([ray_path[0][0]], [ray_path[0][1]], c='r', s=25, label='Start')

        self.ax.legend(loc='best')
        self.ax.grid(True, linestyle='--', alpha=0.3)
        self.ax.set_xlabel("X")
        self.ax.set_ylabel("Y")
        # Auto-zoom based on geometry with margin
        self._set_axes_from_geometry(polygon, margin=10.0)
        self.canvas.draw()

    def _set_axes_from_geometry(self, polygon: Polygon, margin: float = 10.0) -> None:
        xs = [float(v[0]) for v in polygon.vertices]
        ys = [float(v[1]) for v in polygon.vertices]
        min_x, max_x = min(xs), max(xs)
        min_y, max_y = min(ys), max(ys)
        self.ax.set_xlim(min_x - margin, max_x + margin)
        self.ax.set_ylim(min_y - margin, max_y + margin)

    def _boundary_name(self, idx: int) -> str:
        shape = self.shape_var.get()
        if shape == "Rectangle":
            mapping = {
                0: "Bottom",
                1: "Right",
                2: "Top",
                3: "Left",
            }
        else:
            mapping = {
                0: "Bottom",
                1: "Left-slanted",
                2: "Top",
                3: "Right-slanted",
            }
        return mapping.get(idx, "Edge")

    def _update_boundary_mapping(self, polygon: Polygon) -> None:
        shape = self.shape_var.get()
        if shape == "Rectangle":
            text = (
                "Boundary indices:\n"
                "0 = Bottom, 1 = Right, 2 = Top, 3 = Left"
            )
        else:
            text = (
                "Boundary indices:\n"
                "0 = Bottom, 1 = Left-slanted, 2 = Top, 3 = Right-slanted"
            )
        if hasattr(self, "boundary_map_label"):
            self.boundary_map_label.configure(text=text)

    def _toggle_info(self) -> None:
        show = self.show_info_var.get()
        if show:
            self.info_frame.grid()
        else:
            self.info_frame.grid_remove()

    def _boundary_index_map(self, polygon: Polygon) -> list[int]:
        """Map UI boundary order (0..3) to polygon edge indices.

        Rectangle: identity mapping.
        Trapezoid: 0=bottom, 1=left-slanted, 2=top, 3=right-slanted → polygon edges [0,3,2,1].
        """
        shape = self.shape_var.get()
        if shape == "Rectangle":
            return [0, 1, 2, 3]
        # Trapezoid
        return [0, 3, 2, 1]


