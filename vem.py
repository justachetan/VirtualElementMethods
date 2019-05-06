import scipy.io
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable


def mod_wrap(x, a):
    return (x % a)


def square_domain_boundary_condition(points):
    x = points[:, 0]
    y = points[:, 1]

    g = (1 - x) * y * np.sin(np.pi * x)

    return g


def square_domain_rhs(point):
    x = point[0]
    y = point[1]

    f = 15 * np.sin(np.pi * x) * np.sin(np.pi * y)
    return f


def L_domain_rhs(point):
    return 0


def L_domain_boundary_condition(points):
    x = points[:, 0]
    y = points[:, 1]

    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    theta = ((theta >= 0) * theta) + ((theta < 0) * (theta + (2 * np.pi)))
    g = (r**(2 / 3)) * np.sin(2 * (theta - np.pi / 2) / 3)
    return g


def vem(mesh_file, rhs, boundary_condition):

    mesh = scipy.io.loadmat(mesh_file)

    vertices = mesh['vertices']

    elements = np.array([i[0].reshape(-1) - 1 for i in mesh['elements']])

    boundary = mesh['boundary'].T[0] - 1

    n_dofs = vertices.shape[0]
    n_polys = 3

    K = np.zeros((n_dofs, n_dofs))
    F = np.zeros(n_dofs)
    u = np.zeros(n_dofs)

    linear_polynomials = [[0, 0], [1, 0], [0, 1]]

    for el_id in range(elements.shape[0]):

        vert_ids = elements[el_id]

        verts = vertices[vert_ids]

        n_sides = vert_ids.shape[0]

        area_components = verts[
            :, 0] * np.roll(verts[:, 1], -1) - np.roll(verts[:, 0], -1) * verts[:, 1]

        area = 0.5 * np.abs(np.sum(area_components))

        centroid = np.sum((np.roll(verts, -1, axis=0) + verts)
                          * area_components.reshape(-1, 1), axis=0) / (6 * area)

        diameter = np.max(np.linalg.norm(
            verts - np.roll(verts, -1, axis=0), ord=2))

        D = np.zeros((n_sides, n_polys))
        D[:, 0] = 1

        B = np.zeros((n_polys, n_sides))
        B[0, :] = 1 / n_sides

        for vertex_id in range(n_sides):

            vert = verts[vertex_id, :]

            prevv = verts[mod_wrap(vertex_id - 1, n_sides), :]

            nextv = verts[mod_wrap(vertex_id + 1, n_sides), :]

            vertex_normal = np.array(
                [nextv[1] - prevv[1], prevv[0] - nextv[0]])

            for poly_id in range(1, n_polys):  # Looping over non-constant polynomials

                poly_degree = linear_polynomials[poly_id]

                monomial_grad = poly_degree / diameter

                D[vertex_id, poly_id] = np.dot(
                    vert - centroid, poly_degree) / diameter
                B[poly_id, vertex_id] = np.dot(
                    monomial_grad, vertex_normal) / 2

        projector = np.dot(np.linalg.inv(np.dot(B, D)), B)

        stabilising_term = np.dot(
            (np.eye(n_sides) - np.dot(D, projector)).T, (np.eye(n_sides) - np.dot(D, projector)))

        G = np.dot(B, D)
        G[0, :] = 0

        local_stiffness = np.dot(
            np.dot(projector.T, G), projector) + stabilising_term

        # Global indices
        gis = np.array(np.meshgrid(vert_ids, vert_ids)
                       ).T.reshape(-1, 2).tolist()

        lsr = local_stiffness.ravel()
        counter = 0
        for i in range(len(gis)):
            x = gis[i][0]
            y = gis[i][1]
            K[x, y] = K[x, y] + lsr[counter]
            counter += 1

        F[vert_ids] = F[vert_ids] + (rhs(centroid) * (area / n_sides))

    boundary_vals = boundary_condition(vertices[boundary])
    internal_dofs = np.array(
        [i for i in np.arange(n_dofs) if i not in boundary])

    F = F - np.dot(K[:, boundary], boundary_vals)

    num_idof = internal_dofs.shape[0]
    gid_of_idof = np.array(np.meshgrid(
        internal_dofs, internal_dofs)).T.reshape(-1, 2).tolist()
    K_II = np.zeros((num_idof, num_idof)).ravel()
    counter = 0
    for i in range(len(gid_of_idof)):
        x = gid_of_idof[i][0]
        y = gid_of_idof[i][1]

        K_II[counter] = K[x, y]
        counter += 1

    K_II = K_II.reshape(num_idof, num_idof)
    u[internal_dofs] = np.linalg.solve(K_II, F[internal_dofs])
    u[boundary] = boundary_vals

    return u


def plot_solution(mesh_file, u, save=False, plot_name=None):

    mesh = scipy.io.loadmat(mesh_file)

    vertices = mesh['vertices']
    elements = np.array([i[0].reshape(-1) - 1 for i in mesh['elements']])
    boundary = mesh['boundary'].T[0] - 1

    x = vertices[:, 0]
    y = vertices[:, 1]
    v = u

    plt.figure(figsize=(5, 5))
    ax = plt.subplot(111)
    xi = np.linspace(min(x) - 0.01, max(x) + 0.001, 100)
    yi = np.linspace(min(y) - 0.01, max(y) + 0.001, 100)

    zi = griddata((x, y), v, (xi[None, :], yi[:, None]), method='linear')
    for i in range(len(elements)):

        for j in range(len(elements[i])):

            x = [vertices[elements[i][j % len(elements[i])]][0], vertices[
                elements[i][(j + 1) % len(elements[i])]][0]]
            y = [vertices[elements[i][j % len(elements[i])]][1], vertices[
                elements[i][(j + 1) % len(elements[i])]][1]]
            plt.plot(x, y, "k", linewidth=0.5)

    im = plt.pcolormesh(xi, yi, zi)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    ax.set_title("Approximate Solution (u)")
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if save and plot_name is not None:
        plt.savefig(plot_name)
    elif save and plot_name is None:
        plt.savefig("sol.png")

    plt.show()


def main():
    import argparse
    from argparse import RawTextHelpFormatter

    parser = argparse.ArgumentParser(
        description='This command is used for solving 2-D PDEs using Virtual Element Methods\n of the lowest order.', formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        "-d", help="Specifies the shape of the 2D domain.\n Possible values are:\n- s: Square Domain\n- l: L-Shaped Domain", type=str)
    parser.add_argument("-o", help="Path to output file",
                        type=str, default="./sol.npy")
    parser.add_argument("i", help="Path to input mesh", type=str)
    parser.add_argument(
        "--save_plot", help="Flag for saving the plot", action="store_true")
    parser.add_argument("--title", help="Title of plot",
                        default="./plot.png", type=str)

    args = parser.parse_args()

    mesh_file = args.i

    u = None

    if args.d == "s":
        u = vem(mesh_file, square_domain_rhs, square_domain_boundary_condition)
    elif args.d == "l":
        u = vem(mesh_file, L_domain_rhs, L_domain_boundary_condition)
    else:
        raise RuntimeError("Shape of domain not supported!")

    np.save(args.o, u)

    plot_solution(mesh_file, u, args.save_plot, args.title)


if __name__ == '__main__':
    main()
