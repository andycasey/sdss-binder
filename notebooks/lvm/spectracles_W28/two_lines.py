import jax
from equinox import Module
from jaxtyping import Array
from spectracles import (
    AnyParameter,
    Constant,
    ConstrainedParameter,
    Kernel,
    PerSpaxel,
    SpatialDataLVM,
    SpectralSpatialModel,
)
from spectracles.lvm_models.fields import GPField, PositiveGPField
from spectracles.lvm_models.likelihood import ln_likelihood
from spectracles.lvm_models.line_single import EmissionLine


class x2_from_x1_dx(Module):
    x1: AnyParameter
    dx: AnyParameter

    @property
    def val(self):
        return self.x1.val + self.dx.val

    def __call__(self) -> Array:
        return self.x1.val + self.dx.val


class TwoComponentEmissionLine(SpectralSpatialModel):
    # Line components
    line_1: EmissionLine
    line_2: EmissionLine
    # Continuum
    offs: PerSpaxel  # Nuisance offsets per spaxel

    def __init__(
        self,
        n_spaxels: int,
        offsets: AnyParameter,
        line_centre: AnyParameter,
        n_modes: tuple[int, int],
        A_kernel: Kernel,
        v_kernel: Kernel,
        σ_kernel: Kernel,
        σ_lsf: AnyParameter,
        v_bary: AnyParameter,
        v_syst: AnyParameter,
        Δv_syst: ConstrainedParameter,
    ):
        self.offs = Constant(
            const=PerSpaxel(n_spaxels=n_spaxels, spaxel_values=offsets)
        )
        self.line_1 = EmissionLine(
            μ=line_centre,
            A=PositiveGPField(kernel=A_kernel, n_modes=n_modes),
            v=GPField(kernel=v_kernel, n_modes=n_modes),
            vσ=PositiveGPField(kernel=σ_kernel, n_modes=n_modes),
            σ_lsf=PerSpaxel(n_spaxels=n_spaxels, spaxel_values=σ_lsf),
            v_bary=PerSpaxel(n_spaxels=n_spaxels, spaxel_values=v_bary),
            v_syst=v_syst,
        )
        self.line_2 = EmissionLine(
            μ=line_centre,
            A=PositiveGPField(kernel=A_kernel, n_modes=n_modes),
            v=GPField(kernel=v_kernel, n_modes=n_modes),
            vσ=PositiveGPField(kernel=σ_kernel, n_modes=n_modes),
            σ_lsf=PerSpaxel(n_spaxels=n_spaxels, spaxel_values=σ_lsf),
            v_bary=PerSpaxel(n_spaxels=n_spaxels, spaxel_values=v_bary),
            v_syst=x2_from_x1_dx(x1=v_syst, dx=Δv_syst),
        )

    def __call__(self, λ: Array, spatial_data: SpatialDataLVM) -> Array:
        return (
            self.line_1(λ, spatial_data)
            + self.line_2(λ, spatial_data)
            + self.offs(λ, spatial_data)
        )


def neg_ln_posterior(model, λ, xy_data, data, u_data, mask):
    vmapped_model = jax.vmap(model, in_axes=(0, None))
    ln_like = ln_likelihood(vmapped_model, λ, xy_data, data, u_data, mask)
    locked_model = model.get_locked_model()
    ln_prior = (
        locked_model.line_1.A.gp.prior_logpdf()
        + locked_model.line_1.v.gp.prior_logpdf()
        + locked_model.line_1.vσ.gp.prior_logpdf()
        + locked_model.line_2.A.gp.prior_logpdf()
        + locked_model.line_2.v.gp.prior_logpdf()
        + locked_model.line_2.vσ.gp.prior_logpdf()
    )
    return -1 * (ln_like + ln_prior)
