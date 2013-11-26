package My::Build;
use Module::Build;
@ISA = qw(Module::Build);

		sub ACTION_code {
			my $self = shift;
                        $self->SUPER::ACTION_code();
			$self->do_system(qw(make));
		}
		
1;
