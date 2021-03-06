% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RRRRR1V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR1V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR1V1G2A0_coriolisvec_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 03:35:37
% EndTime: 2020-08-07 03:35:59
% DurationCPUTime: 23.07s
% Computational Cost: add. (105537->754), mult. (210942->1477), div. (16155->16), fcn. (163608->36), ass. (0->609)
t3619 = qJ(2,3) + qJ(3,3);
t3565 = cos(qJ(1,3) + t3619);
t3566 = cos(qJ(1,3) - t3619);
t3547 = t3566 + t3565;
t3620 = qJ(2,2) + qJ(3,2);
t3567 = cos(qJ(1,2) + t3620);
t3568 = cos(qJ(1,2) - t3620);
t3548 = t3568 + t3567;
t3621 = qJ(2,1) + qJ(3,1);
t3569 = cos(qJ(1,1) + t3621);
t3570 = cos(qJ(1,1) - t3621);
t3549 = t3570 + t3569;
t4056 = pkin(2) ^ 2;
t3634 = cos(qJ(3,3));
t3610 = t3634 ^ 2;
t3637 = cos(qJ(3,2));
t3613 = t3637 ^ 2;
t3640 = cos(qJ(3,1));
t3616 = t3640 ^ 2;
t3622 = legFrame(3,2);
t3592 = cos(t3622);
t3902 = 0.2e1 * t3592;
t3623 = legFrame(2,2);
t3593 = cos(t3623);
t3901 = 0.2e1 * t3593;
t3624 = legFrame(1,2);
t3594 = cos(t3624);
t3900 = 0.2e1 * t3594;
t3595 = t3634 * pkin(2);
t3647 = 0.1e1 / pkin(3);
t4055 = (t3595 + pkin(3)) * t3647;
t3597 = t3637 * pkin(2);
t4054 = (t3597 + pkin(3)) * t3647;
t3599 = t3640 * pkin(2);
t4053 = (t3599 + pkin(3)) * t3647;
t3635 = cos(qJ(2,3));
t3611 = t3635 ^ 2;
t3625 = sin(qJ(3,3));
t3909 = t3610 - 0.1e1 / 0.2e1;
t3626 = sin(qJ(2,3));
t3942 = t3626 * t3635;
t3517 = t3909 * t3942 + (t3611 - 0.1e1 / 0.2e1) * t3634 * t3625;
t4052 = 0.4e1 * t3517;
t3638 = cos(qJ(2,2));
t3614 = t3638 ^ 2;
t3628 = sin(qJ(3,2));
t3908 = t3613 - 0.1e1 / 0.2e1;
t3629 = sin(qJ(2,2));
t3934 = t3629 * t3638;
t3518 = t3908 * t3934 + (t3614 - 0.1e1 / 0.2e1) * t3637 * t3628;
t4051 = 0.4e1 * t3518;
t3641 = cos(qJ(2,1));
t3617 = t3641 ^ 2;
t3631 = sin(qJ(3,1));
t3907 = t3616 - 0.1e1 / 0.2e1;
t3632 = sin(qJ(2,1));
t3926 = t3632 * t3641;
t3519 = t3907 * t3926 + (t3617 - 0.1e1 / 0.2e1) * t3640 * t3631;
t4050 = 0.4e1 * t3519;
t3589 = sin(t3622);
t4049 = 0.2e1 * t3589;
t3590 = sin(t3623);
t4048 = 0.2e1 * t3590;
t3591 = sin(t3624);
t4047 = 0.2e1 * t3591;
t3596 = t3635 * pkin(2);
t4046 = 0.2e1 * t3596;
t3598 = t3638 * pkin(2);
t4045 = 0.2e1 * t3598;
t3600 = t3641 * pkin(2);
t4044 = 0.2e1 * t3600;
t4043 = -0.2e1 * t3611;
t4042 = -0.2e1 * t3614;
t4041 = -0.2e1 * t3617;
t3646 = pkin(3) ^ 2;
t4040 = -0.2e1 * t3646;
t4039 = MDP(9) * pkin(1);
t4038 = t3547 / 0.2e1;
t4037 = t3548 / 0.2e1;
t4036 = t3549 / 0.2e1;
t3650 = 0.1e1 / pkin(2);
t4035 = t3650 / 0.2e1;
t4034 = t3650 / 0.4e1;
t3583 = 0.2e1 * t3611 - 0.1e1;
t3584 = 0.2e1 * t3614 - 0.1e1;
t3585 = 0.2e1 * t3617 - 0.1e1;
t4033 = MDP(10) * pkin(1);
t4032 = pkin(1) * t3635;
t4031 = pkin(1) * t3638;
t4030 = pkin(1) * t3641;
t4029 = pkin(3) * t3634;
t4028 = pkin(3) * t3637;
t4027 = pkin(3) * t3640;
t3644 = xDP(2);
t4026 = pkin(3) * t3644;
t3645 = xDP(1);
t4025 = pkin(3) * t3645;
t3921 = t3635 * t3645;
t3922 = t3635 * t3644;
t3627 = sin(qJ(1,3));
t3939 = t3627 * t3645;
t3940 = t3627 * t3644;
t3505 = (-t3626 * t3939 - t3922) * t3902 + (t3626 * t3940 - t3921) * t4049;
t3506 = (-t3626 * t3644 + t3627 * t3921) * t3592 - t3589 * (t3626 * t3645 + t3627 * t3922);
t3643 = xDP(3);
t3481 = t3505 * t3625 + 0.2e1 * t3506 * t3634 + t3547 * t3643;
t3601 = 0.1e1 / t3625;
t4005 = t3481 * t3601;
t3796 = t4005 / 0.2e1;
t3475 = t3650 * t3796;
t3577 = pkin(2) + t4029;
t3889 = t3625 * t4025;
t3890 = t3625 * t4026;
t3636 = cos(qJ(1,3));
t3920 = t3636 * t3643;
t3472 = ((-t3577 * t3939 + t3890) * t3592 + (t3577 * t3940 + t3889) * t3589 - t3577 * t3920) * t3635 + t3626 * ((t3577 * t3644 + t3627 * t3889) * t3592 + (t3577 * t3645 - t3627 * t3890) * t3589 + t3625 * pkin(3) * t3920);
t3913 = t3647 * t3650;
t3813 = t3601 * t3913;
t3771 = t3472 * t3813;
t3469 = t3475 + t3771;
t4024 = t3469 * pkin(3);
t3918 = t3638 * t3645;
t3919 = t3638 * t3644;
t3630 = sin(qJ(1,2));
t3931 = t3630 * t3645;
t3932 = t3630 * t3644;
t3507 = (-t3629 * t3931 - t3919) * t3901 + (t3629 * t3932 - t3918) * t4048;
t3508 = (-t3629 * t3644 + t3630 * t3918) * t3593 - t3590 * (t3629 * t3645 + t3630 * t3919);
t3482 = t3507 * t3628 + 0.2e1 * t3508 * t3637 + t3548 * t3643;
t3604 = 0.1e1 / t3628;
t4003 = t3482 * t3604;
t3794 = t4003 / 0.2e1;
t3476 = t3650 * t3794;
t3578 = pkin(2) + t4028;
t3885 = t3628 * t4025;
t3886 = t3628 * t4026;
t3639 = cos(qJ(1,2));
t3917 = t3639 * t3643;
t3473 = ((-t3578 * t3931 + t3886) * t3593 + (t3578 * t3932 + t3885) * t3590 - t3578 * t3917) * t3638 + t3629 * ((t3578 * t3644 + t3630 * t3885) * t3593 + (t3578 * t3645 - t3630 * t3886) * t3590 + t3628 * pkin(3) * t3917);
t3811 = t3604 * t3913;
t3770 = t3473 * t3811;
t3470 = t3476 + t3770;
t4023 = t3470 * pkin(3);
t3915 = t3641 * t3645;
t3916 = t3641 * t3644;
t3633 = sin(qJ(1,1));
t3923 = t3633 * t3645;
t3924 = t3633 * t3644;
t3509 = (-t3632 * t3923 - t3916) * t3900 + (t3632 * t3924 - t3915) * t4047;
t3510 = (-t3632 * t3644 + t3633 * t3915) * t3594 - t3591 * (t3632 * t3645 + t3633 * t3916);
t3483 = t3509 * t3631 + 0.2e1 * t3510 * t3640 + t3549 * t3643;
t3607 = 0.1e1 / t3631;
t4001 = t3483 * t3607;
t3792 = t4001 / 0.2e1;
t3477 = t3650 * t3792;
t3579 = pkin(2) + t4027;
t3881 = t3631 * t4025;
t3882 = t3631 * t4026;
t3642 = cos(qJ(1,1));
t3914 = t3642 * t3643;
t3474 = ((-t3579 * t3923 + t3882) * t3594 + (t3579 * t3924 + t3881) * t3591 - t3579 * t3914) * t3641 + t3632 * ((t3579 * t3644 + t3633 * t3881) * t3594 + (t3579 * t3645 - t3633 * t3882) * t3591 + t3631 * pkin(3) * t3914);
t3809 = t3607 * t3913;
t3769 = t3474 * t3809;
t3471 = t3477 + t3769;
t4022 = t3471 * pkin(3);
t4021 = t3610 * pkin(3);
t4020 = t3613 * pkin(3);
t4019 = t3616 * pkin(3);
t3580 = pkin(1) + t3596;
t3581 = pkin(1) + t3598;
t3582 = pkin(1) + t3600;
t4018 = -0.2e1 * pkin(2) * pkin(3);
t4017 = MDP(8) * t3650;
t4016 = MDP(15) * t3650;
t3895 = t3634 * t4024;
t3451 = t3796 + t3895;
t3460 = -0.2e1 * t4024;
t3520 = -t3627 * t3643 + (-t3589 * t3644 + t3592 * t3645) * t3636;
t3945 = t3625 * t3635;
t3891 = pkin(3) * t3945;
t3526 = t3577 * t3626 + t3891;
t3553 = -pkin(3) + t3595 + 0.2e1 * t4021;
t3612 = t3636 ^ 2;
t3955 = t3601 * t3634;
t3805 = t3481 * t3955;
t3808 = t3625 * t3942;
t3550 = pkin(3) * cos(t3619) + t3580;
t3541 = 0.1e1 / t3550;
t3978 = t3541 * t3636;
t3837 = t3520 * t3978;
t3946 = t3625 * t3626;
t3892 = pkin(3) * t3946;
t3905 = 0.2e1 * t4029;
t3979 = t3541 * t3627;
t4015 = (((0.4e1 * t3469 * t4021 + t3460 + t3805) * t3611 - 0.2e1 * (t3796 + 0.2e1 * t3895) * t3808 - t3460 + t3460 * t3610) * t3612 - 0.2e1 * t3627 * (t3553 * t3942 + ((pkin(2) + t3905) * t3611 - t3577) * t3625) * t3837 - t3805 + t3460 + t3547 * ((-t3451 * t3635 + t3469 * t3892) * t3636 + t3526 * t3520 * t3979)) * t3481;
t3894 = t3637 * t4023;
t3452 = t3794 + t3894;
t3461 = -0.2e1 * t4023;
t3521 = -t3630 * t3643 + (-t3590 * t3644 + t3593 * t3645) * t3639;
t3937 = t3628 * t3638;
t3887 = pkin(3) * t3937;
t3527 = t3578 * t3629 + t3887;
t3554 = -pkin(3) + t3597 + 0.2e1 * t4020;
t3615 = t3639 ^ 2;
t3953 = t3604 * t3637;
t3803 = t3482 * t3953;
t3807 = t3628 * t3934;
t3551 = pkin(3) * cos(t3620) + t3581;
t3543 = 0.1e1 / t3551;
t3975 = t3543 * t3639;
t3835 = t3521 * t3975;
t3938 = t3628 * t3629;
t3888 = pkin(3) * t3938;
t3904 = 0.2e1 * t4028;
t3976 = t3543 * t3630;
t4014 = (((0.4e1 * t3470 * t4020 + t3461 + t3803) * t3614 - 0.2e1 * (t3794 + 0.2e1 * t3894) * t3807 - t3461 + t3461 * t3613) * t3615 - 0.2e1 * t3630 * (t3554 * t3934 + ((pkin(2) + t3904) * t3614 - t3578) * t3628) * t3835 - t3803 + t3461 + t3548 * ((-t3452 * t3638 + t3470 * t3888) * t3639 + t3527 * t3521 * t3976)) * t3482;
t3893 = t3640 * t4022;
t3453 = t3792 + t3893;
t3462 = -0.2e1 * t4022;
t3522 = -t3633 * t3643 + (-t3591 * t3644 + t3594 * t3645) * t3642;
t3929 = t3631 * t3641;
t3883 = pkin(3) * t3929;
t3528 = t3579 * t3632 + t3883;
t3555 = -pkin(3) + t3599 + 0.2e1 * t4019;
t3618 = t3642 ^ 2;
t3951 = t3607 * t3640;
t3801 = t3483 * t3951;
t3806 = t3631 * t3926;
t3552 = pkin(3) * cos(t3621) + t3582;
t3545 = 0.1e1 / t3552;
t3972 = t3545 * t3642;
t3833 = t3522 * t3972;
t3930 = t3631 * t3632;
t3884 = pkin(3) * t3930;
t3903 = 0.2e1 * t4027;
t3973 = t3545 * t3633;
t4013 = (((0.4e1 * t3471 * t4019 + t3462 + t3801) * t3617 - 0.2e1 * (t3792 + 0.2e1 * t3893) * t3806 - t3462 + t3462 * t3616) * t3618 - 0.2e1 * t3633 * (t3555 * t3926 + ((pkin(2) + t3903) * t3617 - t3579) * t3631) * t3833 - t3801 + t3462 + t3549 * ((-t3453 * t3641 + t3471 * t3884) * t3642 + t3528 * t3522 * t3973)) * t3483;
t3586 = sin(t3619);
t3542 = 0.1e1 / t3550 ^ 2;
t3987 = t3520 * t3542;
t3427 = (t3586 * t4024 + (t3586 * t3472 * t3650 + (t3626 / 0.2e1 + (pkin(2) * t3626 + pkin(3) * t3586) * t4035) * t3481) * t3601) * t3987;
t4012 = t3427 * t3541;
t4011 = t3427 * t3601;
t3587 = sin(t3620);
t3544 = 0.1e1 / t3551 ^ 2;
t3985 = t3521 * t3544;
t3428 = (t3587 * t4023 + (t3587 * t3473 * t3650 + (t3629 / 0.2e1 + (pkin(2) * t3629 + pkin(3) * t3587) * t4035) * t3482) * t3604) * t3985;
t4010 = t3428 * t3543;
t4009 = t3428 * t3604;
t3588 = sin(t3621);
t3546 = 0.1e1 / t3552 ^ 2;
t3983 = t3522 * t3546;
t3429 = (t3588 * t4022 + (t3588 * t3474 * t3650 + (t3632 / 0.2e1 + (pkin(2) * t3632 + pkin(3) * t3588) * t4035) * t3483) * t3607) * t3983;
t4008 = t3429 * t3545;
t4007 = t3429 * t3607;
t4006 = t3481 * t3520;
t4004 = t3482 * t3521;
t4002 = t3483 * t3522;
t3535 = -t3634 * t3635 + t3946;
t3943 = t3626 * t3634;
t3538 = t3943 + t3945;
t3967 = t3589 * t3627;
t3493 = t3535 * t3967 - t3538 * t3592;
t4000 = t3493 * t3601;
t3961 = t3592 * t3627;
t3494 = -t3535 * t3961 - t3538 * t3589;
t3999 = t3494 * t3601;
t3536 = -t3637 * t3638 + t3938;
t3935 = t3629 * t3637;
t3539 = t3935 + t3937;
t3965 = t3590 * t3630;
t3495 = t3536 * t3965 - t3539 * t3593;
t3998 = t3495 * t3604;
t3959 = t3593 * t3630;
t3496 = -t3536 * t3959 - t3539 * t3590;
t3997 = t3496 * t3604;
t3537 = -t3640 * t3641 + t3930;
t3927 = t3632 * t3640;
t3540 = t3927 + t3929;
t3963 = t3591 * t3633;
t3497 = t3537 * t3963 - t3540 * t3594;
t3996 = t3497 * t3607;
t3957 = t3594 * t3633;
t3498 = -t3537 * t3957 - t3540 * t3591;
t3995 = t3498 * t3607;
t3556 = t3577 * t3635;
t3529 = t3556 - t3892;
t3499 = t3526 * t3589 - t3529 * t3961;
t3994 = t3499 * t3601;
t3557 = t3578 * t3638;
t3530 = t3557 - t3888;
t3500 = t3527 * t3590 - t3530 * t3959;
t3993 = t3500 * t3604;
t3558 = t3579 * t3641;
t3531 = t3558 - t3884;
t3501 = t3528 * t3591 - t3531 * t3957;
t3992 = t3501 * t3607;
t3502 = t3526 * t3592 + t3529 * t3967;
t3991 = t3502 * t3601;
t3503 = t3527 * t3593 + t3530 * t3965;
t3990 = t3503 * t3604;
t3504 = t3528 * t3594 + t3531 * t3963;
t3989 = t3504 * t3607;
t3514 = t3520 ^ 2;
t3511 = t3514 * t3542;
t3515 = t3521 ^ 2;
t3512 = t3515 * t3544;
t3516 = t3522 ^ 2;
t3513 = t3516 * t3546;
t3988 = t3520 * t3541;
t3986 = t3521 * t3543;
t3984 = t3522 * t3545;
t3982 = t3529 * t3636;
t3981 = t3530 * t3639;
t3980 = t3531 * t3642;
t3977 = t3542 * t3601;
t3974 = t3544 * t3604;
t3971 = t3546 * t3607;
t3970 = t3547 * t3601;
t3969 = t3548 * t3604;
t3968 = t3549 * t3607;
t3966 = t3589 * t3636;
t3964 = t3590 * t3639;
t3962 = t3591 * t3642;
t3960 = t3592 * t3636;
t3958 = t3593 * t3639;
t3956 = t3594 * t3642;
t3602 = 0.1e1 / t3625 ^ 2;
t3651 = 0.1e1 / t4056;
t3954 = t3602 * t3651;
t3605 = 0.1e1 / t3628 ^ 2;
t3952 = t3605 * t3651;
t3608 = 0.1e1 / t3631 ^ 2;
t3950 = t3608 * t3651;
t3949 = t3611 * t3625;
t3948 = t3614 * t3628;
t3947 = t3617 * t3631;
t3944 = t3626 * t3627;
t3941 = t3627 * t3635;
t3936 = t3629 * t3630;
t3933 = t3630 * t3638;
t3928 = t3632 * t3633;
t3925 = t3633 * t3641;
t3463 = t3475 + t3771 / 0.2e1;
t3559 = pkin(1) - 0.2e1 * t3892;
t3880 = pkin(1) * t3946;
t3702 = -pkin(3) + t3880 + t4021;
t3787 = -t3647 * t3651 / 0.2e1;
t3906 = t3646 - t4056;
t3912 = (t3469 * t3646 + (t3463 * t3905 + t3796) * pkin(2)) * t3602 * t3481 * t3787 + ((t3610 * t4040 + t3634 * t4018 + t3906) * t3611 - t3559 * t3556 + pkin(3) * t3702) * t3813 * t3511;
t3464 = t3476 + t3770 / 0.2e1;
t3560 = pkin(1) - 0.2e1 * t3888;
t3879 = pkin(1) * t3938;
t3701 = -pkin(3) + t3879 + t4020;
t3911 = (t3470 * t3646 + (t3464 * t3904 + t3794) * pkin(2)) * t3605 * t3482 * t3787 + ((t3613 * t4040 + t3637 * t4018 + t3906) * t3614 - t3560 * t3557 + pkin(3) * t3701) * t3811 * t3512;
t3465 = t3477 + t3769 / 0.2e1;
t3561 = pkin(1) - 0.2e1 * t3884;
t3878 = pkin(1) * t3930;
t3700 = -pkin(3) + t3878 + t4019;
t3910 = (t3471 * t3646 + (t3465 * t3903 + t3792) * pkin(2)) * t3608 * t3483 * t3787 + ((t3616 * t4040 + t3640 * t4018 + t3906) * t3617 - t3561 * t3558 + pkin(3) * t3700) * t3809 * t3513;
t3899 = 0.2e1 * t3647;
t3898 = pkin(2) * t3946;
t3897 = pkin(2) * t3938;
t3896 = pkin(2) * t3930;
t3850 = t3472 * t3954;
t3448 = t3469 * t3850;
t3795 = t4005 / 0.4e1;
t3729 = (-0.4e1 * t3627 * ((t3795 + t3895) * t3949 + (t3634 * t3795 + t3909 * t4024) * t3942 - t3625 * t3451 / 0.2e1) * t3636 + (0.2e1 * t3612 - 0.2e1) * (t3553 * t3611 + (t3559 * t3634 - t3898) * t3635 - t3702) * t3988 + t3547 * ((t3451 * t3626 + t3469 * t3891) * t3627 - (pkin(1) + t3529) * t3837)) * t3601 * t3650 * t3988;
t3406 = -t3729 / 0.2e1 - t3954 * t4015 / 0.4e1 + t3448;
t3877 = t3406 * t3955;
t3849 = t3473 * t3952;
t3449 = t3470 * t3849;
t3793 = t4003 / 0.4e1;
t3728 = (-0.4e1 * t3630 * ((t3793 + t3894) * t3948 + (t3637 * t3793 + t3908 * t4023) * t3934 - t3628 * t3452 / 0.2e1) * t3639 + (0.2e1 * t3615 - 0.2e1) * (t3554 * t3614 + (t3560 * t3637 - t3897) * t3638 - t3701) * t3986 + t3548 * ((t3452 * t3629 + t3470 * t3887) * t3630 - (pkin(1) + t3530) * t3835)) * t3604 * t3650 * t3986;
t3407 = -t3728 / 0.2e1 - t3952 * t4014 / 0.4e1 + t3449;
t3876 = t3407 * t3953;
t3848 = t3474 * t3950;
t3450 = t3471 * t3848;
t3791 = t4001 / 0.4e1;
t3727 = (-0.4e1 * t3633 * ((t3791 + t3893) * t3947 + (t3640 * t3791 + t3907 * t4022) * t3926 - t3631 * t3453 / 0.2e1) * t3642 + (0.2e1 * t3618 - 0.2e1) * (t3555 * t3617 + (t3561 * t3640 - t3896) * t3641 - t3700) * t3984 + t3549 * ((t3453 * t3632 + t3471 * t3883) * t3633 - (pkin(1) + t3531) * t3833)) * t3607 * t3650 * t3984;
t3408 = -t3727 / 0.2e1 - t3950 * t4013 / 0.4e1 + t3450;
t3875 = t3408 * t3951;
t3874 = t3535 * t4011;
t3873 = t3538 * t4011;
t3872 = t3580 * t4012;
t3871 = t3427 * t3979;
t3870 = t3427 * t3978;
t3869 = t3626 * t4011;
t3868 = t3635 * t4011;
t3867 = t3536 * t4009;
t3866 = t3539 * t4009;
t3865 = t3581 * t4010;
t3864 = t3428 * t3976;
t3863 = t3428 * t3975;
t3862 = t3629 * t4009;
t3861 = t3638 * t4009;
t3860 = t3537 * t4007;
t3859 = t3540 * t4007;
t3858 = t3582 * t4008;
t3857 = t3429 * t3973;
t3856 = t3429 * t3972;
t3855 = t3632 * t4007;
t3854 = t3641 * t4007;
t3853 = t3469 * t3987;
t3852 = t3470 * t3985;
t3851 = t3471 * t3983;
t3478 = t3481 ^ 2;
t3847 = t3478 * t3541 * t3602;
t3479 = t3482 ^ 2;
t3846 = t3479 * t3543 * t3605;
t3480 = t3483 ^ 2;
t3845 = t3480 * t3545 * t3608;
t3844 = t3636 * t4006;
t3843 = t3639 * t4004;
t3842 = t3642 * t4002;
t3841 = t3514 * t3977;
t3840 = t3515 * t3974;
t3839 = t3516 * t3971;
t3832 = t3601 * t3982;
t3831 = t3604 * t3981;
t3830 = t3607 * t3980;
t3829 = t3541 * t3960;
t3828 = t3541 * t3944;
t3827 = t3626 * t3978;
t3826 = t3541 * t3941;
t3825 = t3635 * t3978;
t3824 = t3543 * t3958;
t3823 = t3543 * t3936;
t3822 = t3629 * t3975;
t3821 = t3543 * t3933;
t3820 = t3638 * t3975;
t3819 = t3545 * t3956;
t3818 = t3545 * t3928;
t3817 = t3632 * t3972;
t3816 = t3545 * t3925;
t3815 = t3641 * t3972;
t3403 = -t3729 + 0.2e1 * t3448 + (-t4015 / 0.2e1 - t3469 * t3472 * t4055) * t3954 + t3912;
t3814 = t3403 * t3955;
t3404 = -t3728 + 0.2e1 * t3449 + (-t4014 / 0.2e1 - t3470 * t3473 * t4054) * t3952 + t3911;
t3812 = t3404 * t3953;
t3405 = -t3727 + 0.2e1 * t3450 + (-t4013 / 0.2e1 - t3471 * t3474 * t4053) * t3950 + t3910;
t3810 = t3405 * t3951;
t3804 = t3406 * t3982;
t3802 = t3407 * t3981;
t3800 = t3408 * t3980;
t3799 = t3478 * t4034;
t3798 = t3479 * t4034;
t3797 = t3480 * t4034;
t3790 = t3970 / 0.2e1;
t3789 = t3969 / 0.2e1;
t3788 = t3968 / 0.2e1;
t3786 = t3636 * t3902;
t3785 = t3639 * t3901;
t3784 = t3642 * t3900;
t3783 = t3427 * t3832;
t3603 = t3626 ^ 2;
t3782 = t3603 * t3870;
t3781 = t3427 * t3827;
t3780 = t3427 * t3825;
t3779 = t3428 * t3831;
t3606 = t3629 ^ 2;
t3778 = t3606 * t3863;
t3777 = t3428 * t3822;
t3776 = t3428 * t3820;
t3775 = t3429 * t3830;
t3609 = t3632 ^ 2;
t3774 = t3609 * t3856;
t3773 = t3429 * t3817;
t3772 = t3429 * t3815;
t3768 = t3481 * t3583 * t3987;
t3767 = t3482 * t3584 * t3985;
t3766 = t3483 * t3585 * t3983;
t3765 = t3517 * t3841;
t3523 = -0.4e1 * t3634 * t3808 + (0.4e1 * t3611 - 0.2e1) * t3610 - t3583;
t3764 = t3523 * t3841;
t3763 = t3547 * t3841;
t3762 = (t4046 + pkin(1)) * t3626 * t3511;
t3761 = t3518 * t3840;
t3524 = -0.4e1 * t3637 * t3807 + (0.4e1 * t3614 - 0.2e1) * t3613 - t3584;
t3760 = t3524 * t3840;
t3759 = t3548 * t3840;
t3758 = (t4045 + pkin(1)) * t3629 * t3512;
t3757 = t3519 * t3839;
t3525 = -0.4e1 * t3640 * t3806 + (0.4e1 * t3617 - 0.2e1) * t3616 - t3585;
t3756 = t3525 * t3839;
t3755 = t3549 * t3839;
t3754 = (t4044 + pkin(1)) * t3632 * t3513;
t3753 = t3942 * t3977;
t3752 = t3934 * t3974;
t3751 = t3926 * t3971;
t3750 = t3406 * t3827;
t3749 = t3406 * t3825;
t3748 = t3407 * t3822;
t3747 = t3407 * t3820;
t3746 = t3408 * t3817;
t3745 = t3408 * t3815;
t3744 = -t3427 * t3970 / 0.2e1;
t3743 = -t3428 * t3969 / 0.2e1;
t3742 = -t3429 * t3968 / 0.2e1;
t3741 = -t3847 / 0.4e1;
t3740 = t3847 / 0.4e1;
t3739 = -t3846 / 0.4e1;
t3738 = t3846 / 0.4e1;
t3737 = -t3845 / 0.4e1;
t3736 = t3845 / 0.4e1;
t3735 = -0.2e1 * t3765;
t3734 = -0.2e1 * t3761;
t3733 = -0.2e1 * t3757;
t3732 = -0.2e1 * t3427 * t3826;
t3731 = -0.2e1 * t3428 * t3821;
t3730 = -0.2e1 * t3429 * t3816;
t3487 = (t4043 + 0.1e1) * t3511;
t3488 = (t4042 + 0.1e1) * t3512;
t3489 = (t4041 + 0.1e1) * t3513;
t3726 = t3636 * t3768;
t3725 = t3639 * t3767;
t3724 = t3642 * t3766;
t3723 = t3832 * t3511;
t3722 = t3831 * t3512;
t3721 = t3830 * t3513;
t3720 = t3636 * t3741;
t3719 = t3636 * t3740;
t3718 = t3639 * t3739;
t3717 = t3639 * t3738;
t3716 = t3642 * t3737;
t3715 = t3642 * t3736;
t3400 = -t3448 * t4055 + t3406 + t3912;
t3466 = t3469 ^ 2;
t3714 = t3400 * t3538 - t3466 * t3535;
t3713 = t3400 * t3535 + t3466 * t3538;
t3401 = -t3449 * t4054 + t3407 + t3911;
t3467 = t3470 ^ 2;
t3712 = t3401 * t3539 - t3467 * t3536;
t3711 = t3401 * t3536 + t3467 * t3539;
t3402 = -t3450 * t4053 + t3408 + t3910;
t3468 = t3471 ^ 2;
t3710 = t3402 * t3540 - t3468 * t3537;
t3709 = t3402 * t3537 + t3468 * t3540;
t3708 = -0.2e1 * t3589 * t3780;
t3707 = -0.2e1 * t3590 * t3776;
t3706 = -0.2e1 * t3591 * t3772;
t3705 = t3780 * t3902;
t3704 = t3776 * t3901;
t3703 = t3772 * t3900;
t3699 = pkin(2) * t3463 * t4043 - t3469 * t4032;
t3698 = pkin(2) * t3464 * t4042 - t3470 * t4031;
t3697 = pkin(2) * t3465 * t4041 - t3471 * t4030;
t3696 = t3541 * t3714;
t3695 = t3541 * t3713;
t3694 = t3543 * t3712;
t3693 = t3543 * t3711;
t3692 = t3545 * t3710;
t3691 = t3545 * t3709;
t3690 = t3514 * t4038 + t3627 * t4006;
t3689 = t3515 * t4037 + t3630 * t4004;
t3688 = t3516 * t4036 + t3633 * t4002;
t3687 = (pkin(1) * t3469 + t3463 * t4046) * t3626;
t3686 = (pkin(1) * t3470 + t3464 * t4045) * t3629;
t3685 = (pkin(1) * t3471 + t3465 * t4044) * t3632;
t3684 = t3493 * t3514 + t3589 * t3844;
t3683 = -t3494 * t3514 + t3592 * t3844;
t3682 = t3495 * t3515 + t3590 * t3843;
t3681 = -t3496 * t3515 + t3593 * t3843;
t3680 = t3497 * t3516 + t3591 * t3842;
t3679 = -t3498 * t3516 + t3594 * t3842;
t3678 = -(t3625 * t3687 + (t3796 + t3699) * t3634) * t3987 + t3538 * t3872;
t3677 = -(t3628 * t3686 + (t3794 + t3698) * t3637) * t3985 + t3539 * t3865;
t3676 = -(t3631 * t3685 + (t3792 + t3697) * t3640) * t3983 + t3540 * t3858;
t3675 = ((-t3568 / 0.2e1 - t3567 / 0.2e1) * t3643 + (-t3508 + t3686) * t3637 + (-t3507 / 0.2e1 - t3698) * t3628) * t3985 + t3536 * t3865;
t3674 = ((-t3566 / 0.2e1 - t3565 / 0.2e1) * t3643 + (-t3506 + t3687) * t3634 + (-t3505 / 0.2e1 - t3699) * t3625) * t3987 + t3535 * t3872;
t3673 = ((-t3570 / 0.2e1 - t3569 / 0.2e1) * t3643 + (-t3510 + t3685) * t3640 + (-t3509 / 0.2e1 - t3697) * t3631) * t3983 + t3537 * t3858;
t3672 = 0.2e1 * t3523 * t3853 + t4012 * t4052;
t3671 = t3538 ^ 2 * t4012 + t3853 * t4052;
t3670 = 0.2e1 * t3524 * t3852 + t4010 * t4051;
t3669 = t3539 ^ 2 * t4010 + t3852 * t4051;
t3668 = 0.2e1 * t3525 * t3851 + t4008 * t4050;
t3667 = t3540 ^ 2 * t4008 + t3851 * t4050;
t3666 = t3690 * t3977;
t3665 = t3689 * t3974;
t3664 = t3688 * t3971;
t3663 = t3684 * t3977;
t3662 = t3683 * t3977;
t3661 = t3682 * t3974;
t3660 = t3681 * t3974;
t3659 = t3680 * t3971;
t3658 = t3679 * t3971;
t3657 = 0.2e1 * t3678;
t3656 = 0.2e1 * t3677;
t3655 = 0.2e1 * t3676;
t3654 = 0.2e1 * t3675;
t3653 = 0.2e1 * t3674;
t3652 = 0.2e1 * t3673;
t3648 = 0.1e1 / pkin(3) ^ 2;
t3459 = t3607 * t3797 + (pkin(2) * t3947 + (pkin(1) * t3631 + pkin(2) * t3927) * t3641 + pkin(1) * t3927) * t3513;
t3458 = t3604 * t3798 + (pkin(2) * t3948 + (pkin(1) * t3628 + pkin(2) * t3935) * t3638 + pkin(1) * t3935) * t3512;
t3457 = t3601 * t3799 + (pkin(2) * t3949 + (pkin(1) * t3625 + pkin(2) * t3943) * t3635 + pkin(1) * t3943) * t3511;
t3456 = t3640 * t3608 * t3797 + (t3617 * t3599 + (pkin(1) * t3640 - t3896) * t3641 - t3878) * t3513;
t3455 = t3637 * t3605 * t3798 + (t3614 * t3597 + (pkin(1) * t3637 - t3897) * t3638 - t3879) * t3512;
t3454 = t3634 * t3602 * t3799 + (t3611 * t3595 + (pkin(1) * t3634 - t3898) * t3635 - t3880) * t3511;
t3447 = t3513 * t4030 + ((-t3474 * t3648 - t3483 * t3647) * t3848 - t3489) * pkin(2);
t3446 = t3512 * t4031 + ((-t3473 * t3648 - t3482 * t3647) * t3849 - t3488) * pkin(2);
t3445 = t3511 * t4032 + ((-t3472 * t3648 - t3481 * t3647) * t3850 - t3487) * pkin(2);
t3444 = t3447 * t3640 - t3631 * t3754;
t3443 = t3447 * t3631 + t3640 * t3754;
t3442 = t3446 * t3637 - t3628 * t3758;
t3441 = t3446 * t3628 + t3637 * t3758;
t3440 = t3445 * t3634 - t3625 * t3762;
t3439 = t3445 * t3625 + t3634 * t3762;
t1 = [(t3427 * t3829 + t3428 * t3824 + t3429 * t3819) * MDP(1) + (t3592 * t3782 + t3593 * t3778 + t3594 * t3774 + (t3679 * t3751 + t3681 * t3752 + t3683 * t3753) * t3650) * MDP(4) + (t3626 * t3705 + t3629 * t3704 + t3632 * t3703 + ((t3489 * t3498 + t3594 * t3724) * t3607 + (t3488 * t3496 + t3593 * t3725) * t3604 + (t3487 * t3494 + t3592 * t3726) * t3601) * t3650) * MDP(5) + (-t3592 * t3750 - t3593 * t3748 - t3594 * t3746 + (t3592 * t3635 * t3720 + t3593 * t3638 * t3718 + t3594 * t3641 * t3716) * t3651 + (-t3494 * t3869 - t3496 * t3862 - t3498 * t3855) * t3650) * MDP(6) + (-t3592 * t3749 - t3593 * t3747 - t3594 * t3745 + (t3592 * t3626 * t3719 + t3593 * t3629 * t3717 + t3594 * t3632 * t3715) * t3651 + (-t3494 * t3868 - t3496 * t3861 - t3498 * t3854) * t3650) * MDP(7) + (t3406 * t3999 + t3407 * t3997 + t3408 * t3995) * t4017 + (t3705 + t3704 + t3703 + (-t3626 * t3662 - t3629 * t3660 - t3632 * t3658) * t3650) * t4039 + (-0.2e1 * t3592 * t3781 - 0.2e1 * t3593 * t3777 - 0.2e1 * t3594 * t3773 + (-t3635 * t3662 - t3638 * t3660 - t3641 * t3658) * t3650) * t4033 + (t3667 * t3956 + t3669 * t3958 + t3671 * t3960 + (t3494 * t3735 + t3496 * t3734 + t3498 * t3733 + (-t3499 * t3765 - t3500 * t3761 - t3501 * t3757) * t3899) * t3650) * MDP(11) + (t3668 * t3956 + t3670 * t3958 + t3672 * t3960 + (-t3494 * t3764 - t3496 * t3760 - t3498 * t3756 + (-t3499 * t3764 - t3500 * t3760 - t3501 * t3756) * t3647) * t3650) * MDP(12) + (-t3710 * t3819 - t3712 * t3824 - t3714 * t3829 + (-t3494 * t3873 - t3496 * t3866 - t3498 * t3859 + (-t3499 * t3873 - t3500 * t3866 - t3501 * t3859) * t3647) * t3650) * MDP(13) + (t3709 * t3819 + t3711 * t3824 + t3713 * t3829 + (t3494 * t3874 + t3496 * t3867 + t3498 * t3860 + (t3499 * t3874 + t3500 * t3867 + t3501 * t3860) * t3647) * t3650) * MDP(14) + (t3400 * t3999 + t3401 * t3997 + t3402 * t3995 + (t3400 * t3994 + t3401 * t3993 + t3402 * t3992) * t3647) * t4016 + (t3494 * t3814 + t3496 * t3812 + t3498 * t3810 - t3673 * t3784 - t3675 * t3785 - t3674 * t3786 + (t3439 * t3999 + t3441 * t3997 + t3443 * t3995) * t3650 + (t3499 * t3877 + t3500 * t3876 + t3501 * t3875 + (t3457 * t3994 + t3458 * t3993 + t3459 * t3992) * t3650) * t3647) * MDP(16) + (-t3494 * t3403 - t3496 * t3404 - t3498 * t3405 - t3676 * t3784 - t3677 * t3785 - t3678 * t3786 + (t3440 * t3999 + t3442 * t3997 + t3444 * t3995) * t3650 + (-t3499 * t3406 - t3500 * t3407 - t3501 * t3408 + (t3454 * t3994 + t3455 * t3993 + t3456 * t3992) * t3650) * t3647) * MDP(17); (-t3589 * t3870 - t3590 * t3863 - t3591 * t3856) * MDP(1) + (-t3589 * t3782 - t3590 * t3778 - t3591 * t3774 + (-t3680 * t3751 - t3682 * t3752 - t3684 * t3753) * t3650) * MDP(4) + (t3626 * t3708 + t3629 * t3707 + t3632 * t3706 + ((t3489 * t3497 - t3591 * t3724) * t3607 + (t3488 * t3495 - t3590 * t3725) * t3604 + (t3487 * t3493 - t3589 * t3726) * t3601) * t3650) * MDP(5) + (t3589 * t3750 + t3590 * t3748 + t3591 * t3746 + (t3589 * t3635 * t3719 + t3590 * t3638 * t3717 + t3591 * t3641 * t3715) * t3651 + (-t3493 * t3869 - t3495 * t3862 - t3497 * t3855) * t3650) * MDP(6) + (t3589 * t3749 + t3590 * t3747 + t3591 * t3745 + (t3589 * t3626 * t3720 + t3590 * t3629 * t3718 + t3591 * t3632 * t3716) * t3651 + (-t3493 * t3868 - t3495 * t3861 - t3497 * t3854) * t3650) * MDP(7) + (t3406 * t4000 + t3407 * t3998 + t3408 * t3996) * t4017 + (t3708 + t3707 + t3706 + (t3626 * t3663 + t3629 * t3661 + t3632 * t3659) * t3650) * t4039 + (t3781 * t4049 + t3777 * t4048 + t3773 * t4047 + (t3635 * t3663 + t3638 * t3661 + t3641 * t3659) * t3650) * t4033 + (-t3667 * t3962 - t3669 * t3964 - t3671 * t3966 + (t3493 * t3735 + t3495 * t3734 + t3497 * t3733 + (-t3502 * t3765 - t3503 * t3761 - t3504 * t3757) * t3899) * t3650) * MDP(11) + (-t3668 * t3962 - t3670 * t3964 - t3672 * t3966 + (-t3493 * t3764 - t3495 * t3760 - t3497 * t3756 + (-t3502 * t3764 - t3503 * t3760 - t3504 * t3756) * t3647) * t3650) * MDP(12) + (t3692 * t3962 + t3694 * t3964 + t3696 * t3966 + (-t3493 * t3873 - t3495 * t3866 - t3497 * t3859 + (-t3502 * t3873 - t3503 * t3866 - t3504 * t3859) * t3647) * t3650) * MDP(13) + (-t3691 * t3962 - t3693 * t3964 - t3695 * t3966 + (t3493 * t3874 + t3495 * t3867 + t3497 * t3860 + (t3502 * t3874 + t3503 * t3867 + t3504 * t3860) * t3647) * t3650) * MDP(14) + (t3400 * t4000 + t3401 * t3998 + t3402 * t3996 + (t3400 * t3991 + t3401 * t3990 + t3402 * t3989) * t3647) * t4016 + (t3493 * t3814 + t3495 * t3812 + t3497 * t3810 + t3652 * t3962 + t3654 * t3964 + t3653 * t3966 + (t3439 * t4000 + t3441 * t3998 + t3443 * t3996) * t3650 + (t3502 * t3877 + t3503 * t3876 + t3504 * t3875 + (t3457 * t3991 + t3458 * t3990 + t3459 * t3989) * t3650) * t3647) * MDP(16) + (-t3493 * t3403 - t3495 * t3404 - t3497 * t3405 + t3655 * t3962 + t3656 * t3964 + t3657 * t3966 + (t3440 * t4000 + t3442 * t3998 + t3444 * t3996) * t3650 + (-t3502 * t3406 - t3503 * t3407 - t3504 * t3408 + (t3454 * t3991 + t3455 * t3990 + t3456 * t3989) * t3650) * t3647) * MDP(17); (-t3857 - t3864 - t3871) * MDP(1) + (-t3603 * t3871 - t3606 * t3864 - t3609 * t3857 + (-t3688 * t3751 - t3689 * t3752 - t3690 * t3753) * t3650) * MDP(4) + (t3626 * t3732 + t3629 * t3731 + t3632 * t3730 + ((t3489 * t4036 - t3633 * t3766) * t3607 + (t3488 * t4037 - t3630 * t3767) * t3604 + (t3487 * t4038 - t3627 * t3768) * t3601) * t3650) * MDP(5) + (t3406 * t3828 + t3407 * t3823 + t3408 * t3818 + (t3736 * t3925 + t3738 * t3933 + t3740 * t3941) * t3651 + (t3626 * t3744 + t3629 * t3743 + t3632 * t3742) * t3650) * MDP(6) + (t3406 * t3826 + t3407 * t3821 + t3408 * t3816 + (t3737 * t3928 + t3739 * t3936 + t3741 * t3944) * t3651 + (t3635 * t3744 + t3638 * t3743 + t3641 * t3742) * t3650) * MDP(7) + (t3406 * t3790 + t3407 * t3789 + t3408 * t3788) * t4017 + (t3732 + t3731 + t3730 + (t3626 * t3666 + t3629 * t3665 + t3632 * t3664) * t3650) * t4039 + (0.2e1 * t3427 * t3828 + 0.2e1 * t3428 * t3823 + 0.2e1 * t3429 * t3818 + (t3635 * t3666 + t3638 * t3665 + t3641 * t3664) * t3650) * t4033 + (-t3667 * t3633 - t3669 * t3630 - t3671 * t3627 + (-t3517 * t3763 - t3518 * t3759 - t3519 * t3755 + (t3517 * t3723 + t3518 * t3722 + t3519 * t3721) * t3899) * t3650) * MDP(11) + (-t3668 * t3633 - t3670 * t3630 - t3672 * t3627 + (-t3525 * t3755 / 0.2e1 - t3524 * t3759 / 0.2e1 - t3523 * t3763 / 0.2e1 + (t3523 * t3723 + t3524 * t3722 + t3525 * t3721) * t3647) * t3650) * MDP(12) + (t3633 * t3692 + t3630 * t3694 + t3627 * t3696 + (t3540 * t3742 + t3539 * t3743 + t3538 * t3744 + (t3538 * t3783 + t3539 * t3779 + t3540 * t3775) * t3647) * t3650) * MDP(13) + (-t3633 * t3691 - t3630 * t3693 - t3627 * t3695 + (t3429 * t3537 * t3788 + t3428 * t3536 * t3789 + t3427 * t3535 * t3790 + (-t3535 * t3783 - t3536 * t3779 - t3537 * t3775) * t3647) * t3650) * MDP(14) + (t3402 * t3788 + t3401 * t3789 + t3400 * t3790 + (-t3400 * t3832 - t3401 * t3831 - t3402 * t3830) * t3647) * t4016 + (t3640 * t3405 * t3788 + t3637 * t3404 * t3789 + t3634 * t3403 * t3790 + t3633 * t3652 + t3630 * t3654 + t3627 * t3653 + (t3439 * t3790 + t3441 * t3789 + t3443 * t3788) * t3650 + (-t3804 * t3955 - t3802 * t3953 - t3800 * t3951 + (-t3457 * t3832 - t3458 * t3831 - t3459 * t3830) * t3650) * t3647) * MDP(16) + (-t3549 * t3405 / 0.2e1 - t3548 * t3404 / 0.2e1 - t3547 * t3403 / 0.2e1 + t3633 * t3655 + t3630 * t3656 + t3627 * t3657 + (t3440 * t3790 + t3442 * t3789 + t3444 * t3788) * t3650 + (t3804 + t3802 + t3800 + (-t3454 * t3832 - t3455 * t3831 - t3456 * t3830) * t3650) * t3647) * MDP(17);];
taucX  = t1;
