# import numpy as np
import jax.numpy as np
import scipy.optimize

import jax
from jax.experimental import optimizers
import jax.lax

def minimize_scipy_bfgs(potential_fn, grad_fn, conf):


    N = conf.shape[0]
    D = conf.shape[1]

    def wrap_fun(x):
        x = np.reshape(x, (N, D))
        uu = potential_fn(x)
        gg = grad_fn(x)
        return uu, np.reshape(gg, (-1,))

    res = scipy.optimize.minimize(
        wrap_fun,
        conf,
        jac=True,
        method='BFGS',
        options={'disp': True}
    )

    return np.reshape(res.x, (N, D))

def minimize_jax_adam(potential_fn, grad_fn, conf):

    step_size = 2e-2
    num_epochs = 1000

    # (ytz): Unclear if Adam is actually a sane choice
    opt_init, opt_update, get_params = optimizers.adam(step_size)

    opt_state = opt_init(conf)

    def cond_fun(val):
        idx, _ = val
        return idx < num_epochs

    def body_fun(val):
        epoch, opt_state = val
        x_i = get_params(opt_state)
        dx = grad_fn(x_i)[0]
        opt_state = opt_update(epoch, dx, opt_state)
        return (epoch+1, opt_state)

    init_val = (0, opt_state)

    _, final_state = jax.lax.while_loop(cond_fun, body_fun, init_val)
    return get_params(final_state)

""" Taken from Sam Schoenholz's Jax MD package"""
def minimize_jax_fire(
    potential_fn, grad_fn, conf):
    """Defines FIRE minimization.
    This code implements the "Fast Inertial Relaxation Engine" from [1].
    Args:
    energy_or_force: A function that produces either an energy or a force from
      a set of particle positions specified as an ndarray of shape
      [n, spatial_dimension].
    shift_fn: A function that displaces positions, R, by an amount dR. Both R
      and dR should be ndarrays of shape [n, spatial_dimension].
    quant: Either a quantity.Energy or a quantity.Force specifying whether
      energy_or_force is an energy or force respectively.
    dt_start: The initial step size during minimization as a float.
    dt_max: The maximum step size during minimization as a float.
    n_min: An integer specifying the minimum number of steps moving in the
      correct direction before dt and f_alpha should be updated.
    f_inc: A float specifying the fractional rate by which the step size
      should be increased.
    f_dec: A float specifying the fractional rate by which the step size
      should be decreased.
    alpha_start: A float specifying the initial momentum.
    f_alpha: A float specifying the fractional change in momentum.
    Returns:
    See above.
    [1] Bitzek, Erik, Pekka Koskinen, Franz Gahler, Michael Moseler,
      and Peter Gumbsch. "Structural relaxation made simple."
      Physical review letters 97, no. 17 (2006): 170201.
    """

    dt_start=1e-5
    dt_max=0.1
    n_min=5
    f_inc=1.1
    f_dec=0.5
    alpha_start=0.1
    f_alpha=0.99

    num_epochs = 1000

    def cond_fun(val):
        epoch = val[0]
        return epoch < num_epochs

    def body_fun(val):
        epoch, xi, v, old_f, dt, alpha, n_pos = val
        print("dt", dt, potential_fn(xi))
        xi += dt * v + dt ** 2. * old_f
        f  = -grad_fn(xi)[0]
        # print(f)
        v = v+dt*(old_f+f)/2.
        f_norm = np.sqrt(np.sum(f**2 + 1e-6))
        v_norm = np.sqrt(np.sum(f**2))

        p = np.array(np.dot(np.reshape(f, (-1)), np.reshape(v, (-1))))
        v = v + alpha * (f * v_norm / f_norm - v)

        n_pos = np.where(p >= 0, n_pos + 1.0, 0.)
        dt_choice = np.array([dt * f_inc, dt_max])
        dt = np.where(
            p > 0, np.where(n_pos > n_min, np.min(dt_choice), dt), dt)
        dt = np.where(p < 0, dt * f_dec, dt)
        alpha = np.where(
            p > 0, np.where(n_pos > n_min, alpha * f_alpha, alpha), alpha)
        alpha = np.where(p < 0, alpha_start, alpha)
        v = (p < 0) * np.zeros_like(v) + (p >= 0) * v

        return (epoch+1, xi, v, f, dt, alpha, n_pos)



    init_val = (np.array(0.), conf, np.zeros_like(conf), -grad_fn(conf)[0], np.array(dt_start), np.array(alpha_start), np.array(0.))

    res = jax.lax.while_loop(cond_fun, body_fun, init_val)

    # res = init_val

    # while cond_fun(res):
        # res = body_fun(res)

    return res[1]


