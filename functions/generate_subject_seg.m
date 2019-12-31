function generate_subject_seg(dir_struct,ref_param,mc_param)

% sets up parameters to call generate_bin_vol

% input:
%   dir_struct: structure with fields
%       volume_dir: directory where volume is contained
%       freesurfer_dir: directory where freesurfer can be sourced
%       mr_dir: directory that contains MRI structural scan

%   ref_param: structure with fields
%       has_fiducial: 0 or 1 flag denoting whether volume has fiducial

%   mc_param: structure with fields
%       volume_name: name of volume
%       volume_name_vite: name of volume containing fiducial

%%

if ref_param.has_fiducial
    if ~exist([dir_struct.volume_dir filesep mc_param.volume_name_vite],'file')
        subtext='Creating volume from T1 file';
        print_box('Fiducial-based subject volume does not exist',subtext,70);
        if ~isunix
            error(['Please switch to a Linux workstation so freesurfer can be sourced and MCX can be run%s\n'...
            'Don''t forget to\n'...
            'setenv FREESURFER_HOME /autofs/cluster/freesurfer/centos6_x86_64/stable6_0_0\n'...
            'source $FREESURFER_HOME/SetUpFreeSurfer.csh\n'],'')
        end
        include_fiducial=1;
        generate_bin_vol(dir_struct,mc_param.volume_name,include_fiducial,mc_param.volume_name_vite);
    end
end

if ~exist([dir_struct.volume_dir filesep mc_param.volume_name],'file')
    subtext='Creating volume from T1 file';
    print_box('Subject volume does not exist',subtext,70);
    if ~isunix
        error(['Please switch to a Linux workstation so freesurfer can be sourced and MCX can be run%s\n'...
            'Don''t forget to\n'...
            'setenv FREESURFER_HOME /autofs/cluster/freesurfer/centos6_x86_64/stable6_0_0\n'...
            'source $FREESURFER_HOME/SetUpFreeSurfer.csh\n'],'')
    end
    include_fiducial=0;
    generate_bin_vol(dir_struct,mc_param.volume_name,include_fiducial);
end
    
   