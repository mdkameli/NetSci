function [ param,modified ] = param_init( param )
%PARAM_INIT generates/checks the input parameter struct
%  param = PARAM        returns a default setting
%  param = PARAM(param) returns a valid parameter setting

modified = 0;

if ~nargin,
    param = struct;
    modified = 1;
end;

if ~isfield(param,'function'),
    param.function = @expm;
    disp('Warning: .function not specified, set to @expm.');
    modified = 1;
end;

if ~isfield(param,'thick'),
    param.thick = [];
    disp('Warning: .thick not specified, set to [].');
    modified = 1;
end;


if ~isfield(param,'harmonic_target'),
    param.harmonic_target = inf;
    %disp('Warning: .harmonic_target not specified, set to inf.');
end;

if ~isfield(param,'inner_product'),
    param.inner_product = @inner_product;
    disp('Warning: .inner_product not specified, set to @inner_product.');
    modified = 1;
end;

if ~isfield(param,'V_full'),
    param.V_full = false;
    disp('Warning: .V_full not specified, set to false.');
    modified = 1;
end;

if ~isfield(param,'H_full'),
    param.H_full = true;
    disp('Warning: .H_full not specified, set to true.');
    modified = 1;
end;

if ~param.H_full && ~isstruct(param.function),
    param.H_full = true;
    disp('Warning: .H_full == false only available if .function is parfrac! Set to true.');
    modified = 1;
end;

if param.H_full && isstruct(param.function),
    disp('Warning: You may set .H_full = 0 for faster computation.');
    modified = 1;
end;

if ~isfield(param,'restart_length'),
    param.restart_length = 10;
    disp('Warning: .restart_length not specified, set to 10.');
    modified = 1;
end;

if ~isfield(param,'max_restarts'),
    param.max_restarts = 10;
    disp('Warning: .max_restarts not specified, set to 10.');
    modified = 1;
end;

if ~isfield(param,'hermitian'),
    param.hermitian = false;
    disp('Warning: .hermitian not specified, set to false.');
    modified = 1;
end;

if ~param.hermitian,
    % param.truncation_length = inf;
else
    param.truncation_length = 2;
    param.reorth_number = 0;
end;

if ~isfield(param,'truncation_length'),
    param.truncation_length = inf;
end;

if ~isfield(param,'reorth_number'),
    param.reorth_number = 0;
    disp('Warning: .reorth_number not specified, set to 0.');
    modified = 1;
end;


if ~isfield(param,'exact'),
    param.exact = [];
    disp('Warning: .exact not specified, no error available.');
    modified = 1;
end;

if ~isfield(param,'bound'),
    param.bound = false;
    disp('Warning: .bound not specified, set to false.');
    modified = 1;
end;

if param.bound && ~param.hermitian,
    disp('Warning: Bound only valid in Hermitian case!');
end;

if ~isfield(param,'stopping_accuracy'),
    param.stopping_accuracy = 1e-12;
    disp('Warning: .stopping_accuracy not specified, set to 1e-12.');
    modified = 1;
end;

if ~isfield(param,'min_decay'),
    param.min_decay = 0.95;
    disp('Warning: .min_decay not specified, set to 0.95.');
    modified = 1;
end;

if ~isfield(param,'waitbar'),
    param.waitbar = true;
    disp('Warning: .waitbar not specified, set to true.');
end;
