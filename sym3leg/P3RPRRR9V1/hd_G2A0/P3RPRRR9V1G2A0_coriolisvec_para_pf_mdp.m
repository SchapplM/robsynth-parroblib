% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3RPRRR9V1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G2A0_convert_par2_MPV_fixb.m

% Output:
% taucX [3x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G2A0_coriolisvec_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:52:57
% EndTime: 2020-08-06 18:53:14
% DurationCPUTime: 17.64s
% Computational Cost: add. (69909->594), mult. (111920->1194), div. (6438->21), fcn. (82854->36), ass. (0->504)
t3416 = sin(pkin(7));
t3417 = cos(pkin(7));
t3820 = MDP(4) * t3417 - MDP(5) * t3416;
t3812 = 2 * pkin(3);
t3421 = legFrame(3,2);
t3394 = sin(t3421);
t3397 = cos(t3421);
t3438 = xDP(2);
t3439 = xDP(1);
t3333 = t3394 * t3439 + t3397 * t3438;
t3819 = 0.2e1 * t3333;
t3422 = legFrame(2,2);
t3395 = sin(t3422);
t3398 = cos(t3422);
t3334 = t3395 * t3439 + t3398 * t3438;
t3818 = 0.2e1 * t3334;
t3423 = legFrame(1,2);
t3396 = sin(t3423);
t3399 = cos(t3423);
t3335 = t3396 * t3439 + t3399 * t3438;
t3817 = 0.2e1 * t3335;
t3382 = t3417 * pkin(2);
t3360 = t3382 + pkin(1);
t3411 = pkin(7) + qJ(3,1);
t3381 = cos(t3411);
t3795 = pkin(3) * t3381;
t3350 = t3360 + t3795;
t3410 = pkin(7) + qJ(3,2);
t3380 = cos(t3410);
t3796 = pkin(3) * t3380;
t3349 = t3360 + t3796;
t3409 = pkin(7) + qJ(3,3);
t3379 = cos(t3409);
t3797 = pkin(3) * t3379;
t3348 = t3360 + t3797;
t3425 = sin(qJ(1,3));
t3431 = cos(qJ(1,3));
t3437 = xDP(3);
t3511 = t3394 * t3438 - t3397 * t3439;
t3315 = -t3511 * t3425 + t3431 * t3437;
t3376 = sin(t3409);
t3692 = t3376 * t3333;
t3291 = t3315 * t3379 + t3692;
t3364 = 0.1e1 / t3379;
t3457 = t3379 ^ 2;
t3365 = 0.1e1 / t3457;
t3418 = pkin(5) + qJ(2,3);
t3405 = pkin(6) + t3418;
t3385 = 0.1e1 / t3405;
t3725 = t3333 * t3385;
t3580 = t3365 * t3725;
t3451 = 1 / pkin(3);
t3723 = t3333 * t3451;
t3807 = 0.2e1 * t3409;
t3279 = (-t3291 * sin(t3807) / 0.2e1 + t3348 * t3723) * t3364 * t3580;
t3386 = 0.1e1 / t3405 ^ 2;
t3387 = t3385 * t3386;
t3408 = t3417 ^ 2;
t3430 = cos(qJ(3,3));
t3412 = t3430 ^ 2;
t3424 = sin(qJ(3,3));
t3403 = pkin(1) * t3437;
t3670 = (-t3511 * pkin(1) + t3405 * t3437) * t3425 + (t3511 * t3405 + t3403) * t3431;
t3676 = t3430 * t3424;
t3679 = t3416 * t3424;
t3383 = pkin(1) * t3416;
t3716 = (-pkin(3) * t3424 + t3383) * t3430;
t3724 = t3333 * t3424;
t3737 = t3315 * t3412;
t3264 = ((t3315 * t3430 + t3724) * pkin(2) + (t3676 * t3819 - t3315 + 0.2e1 * t3737) * pkin(3)) * t3408 + (t3670 * t3430 + pkin(1) * t3724 + ((-t3315 * t3424 + t3333 * t3430) * pkin(2) + (-0.2e1 * t3315 * t3676 + t3412 * t3819 - t3333) * pkin(3)) * t3416) * t3417 - pkin(3) * t3737 + t3333 * t3716 - t3670 * t3679 + pkin(3) * t3315;
t3343 = t3417 * t3430 - t3679;
t3337 = 0.1e1 / t3343;
t3767 = t3264 * t3337;
t3619 = t3387 * t3767;
t3552 = t3348 * t3619;
t3710 = t3364 * t3385;
t3611 = t3291 * t3710;
t3286 = pkin(1) * t3611;
t3444 = 0.2e1 * pkin(7);
t3447 = (qJ(2,3) ^ 2);
t3450 = pkin(3) ^ 2;
t3454 = pkin(2) ^ 2;
t3455 = pkin(1) ^ 2;
t3481 = -t3454 * cos(t3444) - (2 * pkin(6) ^ 2) - t3450 - t3454 - 0.2e1 * t3455 + ((-4 * pkin(6) - 2 * pkin(5)) * pkin(5));
t3618 = t3385 * t3767;
t3656 = pkin(2) * t3812;
t3749 = t3291 * t3385;
t3801 = -4 * pkin(5) - 4 * pkin(6);
t3813 = 0.2e1 * pkin(1);
t3646 = (t3618 * t3813 + 0.4e1 * (-t3382 - t3797) * (t3286 - t3618 / 0.2e1) + (0.2e1 * t3405 * t3692 + ((qJ(2,3) * t3801) - t3450 * cos(t3807) - (2 * t3447) + (-cos(qJ(3,3) + t3444) - t3430) * t3656 + t3481) * t3749) * t3364) * t3386 / 0.2e1;
t3751 = t3291 * t3364;
t3474 = -(t3552 + t3646) * t3751 + t3279;
t3427 = sin(qJ(1,2));
t3433 = cos(qJ(1,2));
t3510 = t3395 * t3438 - t3398 * t3439;
t3316 = -t3510 * t3427 + t3433 * t3437;
t3377 = sin(t3410);
t3691 = t3377 * t3334;
t3292 = t3316 * t3380 + t3691;
t3368 = 0.1e1 / t3380;
t3460 = t3380 ^ 2;
t3369 = 0.1e1 / t3460;
t3419 = pkin(5) + qJ(2,2);
t3406 = pkin(6) + t3419;
t3388 = 1 / t3406;
t3722 = t3334 * t3388;
t3578 = t3369 * t3722;
t3720 = t3334 * t3451;
t3806 = 0.2e1 * t3410;
t3280 = (-t3292 * sin(t3806) / 0.2e1 + t3349 * t3720) * t3368 * t3578;
t3389 = 1 / t3406 ^ 2;
t3390 = t3388 * t3389;
t3432 = cos(qJ(3,2));
t3413 = t3432 ^ 2;
t3426 = sin(qJ(3,2));
t3669 = (-t3510 * pkin(1) + t3406 * t3437) * t3427 + (t3510 * t3406 + t3403) * t3433;
t3675 = t3432 * t3426;
t3678 = t3416 * t3426;
t3715 = (-pkin(3) * t3426 + t3383) * t3432;
t3721 = t3334 * t3426;
t3736 = t3316 * t3413;
t3265 = ((t3316 * t3432 + t3721) * pkin(2) + (t3675 * t3818 - t3316 + 0.2e1 * t3736) * pkin(3)) * t3408 + (t3669 * t3432 + pkin(1) * t3721 + ((-t3316 * t3426 + t3334 * t3432) * pkin(2) + (-0.2e1 * t3316 * t3675 + t3413 * t3818 - t3334) * pkin(3)) * t3416) * t3417 - pkin(3) * t3736 + t3334 * t3715 - t3669 * t3678 + pkin(3) * t3316;
t3344 = t3417 * t3432 - t3678;
t3338 = 0.1e1 / t3344;
t3766 = t3265 * t3338;
t3617 = t3390 * t3766;
t3550 = t3349 * t3617;
t3704 = t3368 * t3388;
t3606 = t3292 * t3704;
t3287 = pkin(1) * t3606;
t3448 = qJ(2,2) ^ 2;
t3616 = t3388 * t3766;
t3744 = t3292 * t3388;
t3645 = (t3616 * t3813 + 0.4e1 * (-t3382 - t3796) * (t3287 - t3616 / 0.2e1) + (0.2e1 * t3406 * t3691 + ((qJ(2,2) * t3801) - t3450 * cos(t3806) - (2 * t3448) + (-cos(t3444 + qJ(3,2)) - t3432) * t3656 + t3481) * t3744) * t3368) * t3389 / 0.2e1;
t3746 = t3292 * t3368;
t3473 = -(t3550 + t3645) * t3746 + t3280;
t3429 = sin(qJ(1,1));
t3435 = cos(qJ(1,1));
t3509 = t3396 * t3438 - t3399 * t3439;
t3317 = -t3509 * t3429 + t3435 * t3437;
t3378 = sin(t3411);
t3690 = t3378 * t3335;
t3293 = t3317 * t3381 + t3690;
t3372 = 0.1e1 / t3381;
t3463 = t3381 ^ 2;
t3373 = 0.1e1 / t3463;
t3420 = pkin(5) + qJ(2,1);
t3407 = pkin(6) + t3420;
t3391 = 1 / t3407;
t3719 = t3335 * t3391;
t3576 = t3373 * t3719;
t3717 = t3335 * t3451;
t3805 = 0.2e1 * t3411;
t3281 = (-t3293 * sin(t3805) / 0.2e1 + t3350 * t3717) * t3372 * t3576;
t3392 = 1 / t3407 ^ 2;
t3393 = t3391 * t3392;
t3434 = cos(qJ(3,1));
t3414 = t3434 ^ 2;
t3428 = sin(qJ(3,1));
t3668 = (-t3509 * pkin(1) + t3407 * t3437) * t3429 + (t3509 * t3407 + t3403) * t3435;
t3674 = t3434 * t3428;
t3677 = t3416 * t3428;
t3714 = (-pkin(3) * t3428 + t3383) * t3434;
t3718 = t3335 * t3428;
t3735 = t3317 * t3414;
t3266 = ((t3317 * t3434 + t3718) * pkin(2) + (t3674 * t3817 - t3317 + 0.2e1 * t3735) * pkin(3)) * t3408 + (t3668 * t3434 + pkin(1) * t3718 + ((-t3317 * t3428 + t3335 * t3434) * pkin(2) + (-0.2e1 * t3317 * t3674 + t3414 * t3817 - t3335) * pkin(3)) * t3416) * t3417 - pkin(3) * t3735 + t3335 * t3714 - t3668 * t3677 + pkin(3) * t3317;
t3342 = t3417 * t3434 - t3677;
t3336 = 0.1e1 / t3342;
t3765 = t3266 * t3336;
t3615 = t3393 * t3765;
t3548 = t3350 * t3615;
t3698 = t3372 * t3391;
t3601 = t3293 * t3698;
t3285 = pkin(1) * t3601;
t3449 = qJ(2,1) ^ 2;
t3614 = t3391 * t3765;
t3739 = t3293 * t3391;
t3644 = (t3614 * t3813 + 0.4e1 * (-t3382 - t3795) * (t3285 - t3614 / 0.2e1) + (0.2e1 * t3407 * t3690 + ((qJ(2,1) * t3801) - t3450 * cos(t3805) - (2 * t3449) + (-cos(qJ(3,1) + t3444) - t3434) * t3656 + t3481) * t3739) * t3372) * t3392 / 0.2e1;
t3741 = t3293 * t3372;
t3472 = -(t3548 + t3644) * t3741 + t3281;
t3808 = 0.2e1 * t3408;
t3804 = -0.2e1 * t3424;
t3803 = -0.2e1 * t3426;
t3802 = -0.2e1 * t3428;
t3800 = 0.1e1 - 0.2e1 * t3408;
t3799 = 0.4e1 * t3408 - 0.2e1;
t3798 = pkin(2) * t3408;
t3374 = t3372 * t3373;
t3549 = t3293 * t3615;
t3697 = t3372 * t3417;
t3650 = pkin(2) * t3697;
t3332 = t3335 ^ 2;
t3726 = t3332 * t3391;
t3738 = t3293 * t3392;
t3255 = t3451 * t3374 * t3726 + (-(-t3285 + (t3765 + (-pkin(3) - t3650) * t3293) * t3391) * t3738 - t3549) * t3372;
t3794 = t3255 * pkin(1);
t3366 = t3364 * t3365;
t3553 = t3291 * t3619;
t3709 = t3364 * t3417;
t3652 = pkin(2) * t3709;
t3330 = t3333 ^ 2;
t3728 = t3330 * t3385;
t3748 = t3291 * t3386;
t3256 = t3451 * t3366 * t3728 + (-(-t3286 + (t3767 + (-pkin(3) - t3652) * t3291) * t3385) * t3748 - t3553) * t3364;
t3793 = t3256 * pkin(1);
t3370 = t3368 * t3369;
t3551 = t3292 * t3617;
t3703 = t3368 * t3417;
t3651 = pkin(2) * t3703;
t3331 = t3334 ^ 2;
t3727 = t3331 * t3388;
t3743 = t3292 * t3389;
t3257 = t3451 * t3370 * t3727 + (-(-t3287 + (t3766 + (-pkin(3) - t3651) * t3292) * t3388) * t3743 - t3551) * t3368;
t3792 = t3257 * pkin(1);
t3791 = t3412 * pkin(3);
t3790 = t3413 * pkin(3);
t3789 = t3414 * pkin(3);
t3788 = t3430 * pkin(2);
t3787 = t3432 * pkin(2);
t3786 = t3434 * pkin(2);
t3783 = MDP(12) * t3451 / t3450;
t3244 = t3472 - t3794;
t3782 = t3244 * t3391;
t3246 = t3474 - t3793;
t3781 = t3246 * t3385;
t3248 = t3473 - t3792;
t3780 = t3248 * t3388;
t3314 = t3350 * t3435 + t3407 * t3429;
t3779 = t3255 * t3314;
t3347 = t3416 * t3434 + t3417 * t3428;
t3778 = t3255 * t3347;
t3777 = t3255 * t3372;
t3312 = t3348 * t3431 + t3405 * t3425;
t3776 = t3256 * t3312;
t3345 = t3416 * t3430 + t3417 * t3424;
t3775 = t3256 * t3345;
t3774 = t3256 * t3364;
t3313 = t3349 * t3433 + t3406 * t3427;
t3773 = t3257 * t3313;
t3346 = t3416 * t3432 + t3417 * t3426;
t3772 = t3257 * t3346;
t3771 = t3257 * t3368;
t3770 = (t3286 + (t3291 * t3652 - 0.2e1 * t3767) * t3385) * t3385;
t3769 = (t3287 + (t3292 * t3651 - 0.2e1 * t3766) * t3388) * t3388;
t3768 = (t3285 + (t3293 * t3650 - 0.2e1 * t3765) * t3391) * t3391;
t3558 = pkin(1) * t3425 - t3405 * t3431;
t3297 = t3558 * t3679 + (t3412 - 0.1e1) * t3425 * pkin(3);
t3318 = pkin(1) * t3424 + (-pkin(3) + t3788 + 0.2e1 * t3791) * t3416;
t3445 = -pkin(3) / 0.2e1;
t3351 = t3791 + t3788 / 0.2e1 + t3445;
t3571 = t3425 * t3679;
t3499 = pkin(2) * t3571 + (t3571 * t3812 - t3558) * t3430;
t3686 = t3394 * t3425;
t3446 = pkin(2) / 0.2e1;
t3713 = (pkin(3) * t3430 + t3446) * t3424;
t3273 = (-t3351 * t3686 + t3397 * t3713) * t3808 + (t3397 * t3318 + t3499 * t3394) * t3417 + t3297 * t3394 + t3397 * t3716;
t3764 = t3273 * t3256;
t3683 = t3397 * t3425;
t3274 = (t3351 * t3683 + t3394 * t3713) * t3808 + (t3394 * t3318 - t3499 * t3397) * t3417 - t3297 * t3397 + t3394 * t3716;
t3763 = t3274 * t3256;
t3557 = pkin(1) * t3427 - t3406 * t3433;
t3298 = t3557 * t3678 + (t3413 - 0.1e1) * t3427 * pkin(3);
t3319 = pkin(1) * t3426 + (-pkin(3) + t3787 + 0.2e1 * t3790) * t3416;
t3352 = t3790 + t3787 / 0.2e1 + t3445;
t3570 = t3427 * t3678;
t3498 = pkin(2) * t3570 + (t3570 * t3812 - t3557) * t3432;
t3685 = t3395 * t3427;
t3712 = (pkin(3) * t3432 + t3446) * t3426;
t3275 = (-t3352 * t3685 + t3398 * t3712) * t3808 + (t3398 * t3319 + t3498 * t3395) * t3417 + t3298 * t3395 + t3398 * t3715;
t3762 = t3275 * t3257;
t3682 = t3398 * t3427;
t3276 = (t3352 * t3682 + t3395 * t3712) * t3808 + (t3395 * t3319 - t3498 * t3398) * t3417 - t3298 * t3398 + t3395 * t3715;
t3761 = t3276 * t3257;
t3559 = t3429 * pkin(1) - t3407 * t3435;
t3299 = t3559 * t3677 + (t3414 - 0.1e1) * t3429 * pkin(3);
t3320 = pkin(1) * t3428 + (-pkin(3) + t3786 + 0.2e1 * t3789) * t3416;
t3353 = t3789 + t3786 / 0.2e1 + t3445;
t3569 = t3429 * t3677;
t3497 = pkin(2) * t3569 + (t3569 * t3812 - t3559) * t3434;
t3684 = t3396 * t3429;
t3711 = (pkin(3) * t3434 + t3446) * t3428;
t3277 = (-t3353 * t3684 + t3399 * t3711) * t3808 + (t3399 * t3320 + t3497 * t3396) * t3417 + t3299 * t3396 + t3399 * t3714;
t3760 = t3277 * t3255;
t3681 = t3399 * t3429;
t3278 = (t3353 * t3681 + t3396 * t3711) * t3808 + (t3396 * t3320 - t3497 * t3399) * t3417 - t3299 * t3399 + t3396 * t3714;
t3759 = t3278 * t3255;
t3288 = t3291 ^ 2;
t3758 = t3288 * t3394;
t3757 = t3288 * t3397;
t3289 = t3292 ^ 2;
t3756 = t3289 * t3395;
t3755 = t3289 * t3398;
t3290 = t3293 ^ 2;
t3754 = t3290 * t3396;
t3753 = t3290 * t3399;
t3752 = t3291 * t3333;
t3750 = t3291 * t3365;
t3747 = t3292 * t3334;
t3745 = t3292 * t3369;
t3742 = t3293 * t3335;
t3740 = t3293 * t3373;
t3321 = t3376 * t3397 - t3379 * t3686;
t3734 = t3321 * t3364;
t3322 = t3376 * t3394 + t3379 * t3683;
t3733 = t3322 * t3364;
t3323 = t3377 * t3398 - t3380 * t3685;
t3732 = t3323 * t3368;
t3324 = t3377 * t3395 + t3380 * t3682;
t3731 = t3324 * t3368;
t3325 = t3378 * t3399 - t3381 * t3684;
t3730 = t3325 * t3372;
t3326 = t3378 * t3396 + t3381 * t3681;
t3729 = t3326 * t3372;
t3708 = t3365 * t3387;
t3707 = t3366 * t3376;
t3706 = t3366 * t3386;
t3705 = 0.1e1 / t3457 ^ 2 * t3376;
t3702 = t3369 * t3390;
t3701 = t3370 * t3377;
t3700 = t3370 * t3389;
t3699 = 0.1e1 / t3460 ^ 2 * t3377;
t3696 = t3373 * t3393;
t3695 = t3374 * t3378;
t3694 = t3374 * t3392;
t3693 = 0.1e1 / t3463 ^ 2 * t3378;
t3689 = t3385 * t3431;
t3688 = t3388 * t3433;
t3687 = t3391 * t3435;
t3680 = t3416 * t3417;
t3579 = t3418 * t3723;
t3532 = t3579 / 0.2e1;
t3543 = t3424 * t3611;
t3673 = t3364 * t3430 * t3532 + pkin(1) * t3543;
t3577 = t3419 * t3720;
t3531 = t3577 / 0.2e1;
t3541 = t3426 * t3606;
t3672 = t3368 * t3432 * t3531 + pkin(1) * t3541;
t3575 = t3420 * t3717;
t3530 = t3575 / 0.2e1;
t3539 = t3428 * t3601;
t3671 = t3372 * t3434 * t3530 + pkin(1) * t3539;
t3667 = t3408 - 0.1e1 / 0.2e1;
t3663 = 0.2e1 * t3264 * t3291;
t3662 = 0.2e1 * t3265 * t3292;
t3661 = 0.2e1 * t3266 * t3293;
t3660 = t3386 * t3819;
t3659 = t3389 * t3818;
t3658 = t3392 * t3817;
t3657 = -0.4e1 * t3680;
t3655 = t3255 * t3798;
t3654 = t3256 * t3798;
t3653 = t3257 * t3798;
t3649 = qJ(2,1) * t3696;
t3648 = qJ(2,2) * t3702;
t3647 = qJ(2,3) * t3708;
t3643 = t3336 * t3778;
t3642 = t3255 * t3698;
t3641 = t3396 * t3777;
t3640 = t3399 * t3777;
t3639 = t3420 * t3777;
t3638 = t3255 * t3687;
t3637 = t3337 * t3775;
t3636 = t3256 * t3710;
t3635 = t3394 * t3774;
t3634 = t3397 * t3774;
t3633 = t3418 * t3774;
t3632 = t3256 * t3689;
t3631 = t3338 * t3772;
t3630 = t3257 * t3704;
t3629 = t3395 * t3771;
t3628 = t3398 * t3771;
t3627 = t3419 * t3771;
t3626 = t3257 * t3688;
t3625 = t3394 * t3770;
t3624 = t3397 * t3770;
t3623 = t3395 * t3769;
t3622 = t3398 * t3769;
t3621 = t3396 * t3768;
t3620 = t3399 * t3768;
t3613 = t3321 * t3752;
t3612 = t3322 * t3752;
t3610 = t3365 * t3748;
t3609 = t3430 * t3749;
t3608 = t3323 * t3747;
t3607 = t3324 * t3747;
t3605 = t3369 * t3743;
t3604 = t3432 * t3744;
t3603 = t3325 * t3742;
t3602 = t3326 * t3742;
t3600 = t3373 * t3738;
t3599 = t3434 * t3739;
t3303 = (t3412 - 0.1e1 / 0.2e1) * t3680 + t3667 * t3676;
t3598 = t3303 * t3706;
t3304 = (t3413 - 0.1e1 / 0.2e1) * t3680 + t3667 * t3675;
t3597 = t3304 * t3700;
t3305 = (t3414 - 0.1e1 / 0.2e1) * t3680 + t3667 * t3674;
t3596 = t3305 * t3694;
t3309 = t3799 * t3412 + t3657 * t3676 + t3800;
t3595 = t3309 * t3706;
t3310 = t3799 * t3413 + t3657 * t3675 + t3800;
t3594 = t3310 * t3700;
t3311 = t3799 * t3414 + t3657 * t3674 + t3800;
t3593 = t3311 * t3694;
t3592 = t3321 * t3710;
t3591 = t3322 * t3710;
t3590 = t3323 * t3704;
t3589 = t3324 * t3704;
t3588 = t3325 * t3698;
t3587 = t3326 * t3698;
t3586 = t3330 * t3705;
t3585 = t3330 * t3689;
t3584 = t3331 * t3699;
t3583 = t3331 * t3688;
t3582 = t3332 * t3693;
t3581 = t3332 * t3687;
t3574 = t3364 * t3689;
t3573 = t3368 * t3688;
t3572 = t3372 * t3687;
t3568 = 0.2e1 * t3638;
t3567 = 0.2e1 * t3632;
t3566 = 0.2e1 * t3626;
t3523 = pkin(2) * t3364 * t3609;
t3270 = -t3408 * t3523 + (t3424 * t3532 + (-pkin(1) * t3430 + pkin(2) * t3679) * t3749) * t3709 + t3673 * t3416;
t3565 = 0.2e1 * t3270 * t3725;
t3522 = pkin(2) * t3368 * t3604;
t3271 = -t3408 * t3522 + (t3426 * t3531 + (-pkin(1) * t3432 + pkin(2) * t3678) * t3744) * t3703 + t3672 * t3416;
t3564 = 0.2e1 * t3271 * t3722;
t3521 = pkin(2) * t3372 * t3599;
t3272 = -t3408 * t3521 + (t3428 * t3530 + (-pkin(1) * t3434 + pkin(2) * t3677) * t3739) * t3697 + t3671 * t3416;
t3563 = 0.2e1 * t3272 * t3719;
t3562 = t3337 * t3660;
t3561 = t3338 * t3659;
t3560 = t3336 * t3658;
t3341 = t3347 ^ 2;
t3556 = t3341 * t3642;
t3339 = t3345 ^ 2;
t3555 = t3339 * t3636;
t3340 = t3346 ^ 2;
t3554 = t3340 * t3630;
t3547 = t3288 * t3337 * t3708;
t3546 = t3289 * t3338 * t3702;
t3545 = t3290 * t3336 * t3696;
t3544 = t3312 * t3610;
t3542 = t3313 * t3605;
t3540 = t3314 * t3600;
t3538 = t3394 * t3633;
t3537 = t3395 * t3627;
t3536 = t3396 * t3639;
t3535 = t3397 * t3633;
t3534 = t3398 * t3627;
t3533 = t3399 * t3639;
t3529 = 0.4e1 * t3305 * t3642;
t3528 = 0.4e1 * t3303 * t3636;
t3527 = 0.4e1 * t3304 * t3630;
t3267 = t3543 * t3798 + (t3416 * t3523 + t3673) * t3417 + (pkin(1) * t3609 - t3424 * t3579 / 0.2e1) * t3364 * t3416;
t3526 = -0.2e1 * t3267 * t3580;
t3268 = t3541 * t3798 + (t3416 * t3522 + t3672) * t3417 + (pkin(1) * t3604 - t3426 * t3577 / 0.2e1) * t3368 * t3416;
t3525 = -0.2e1 * t3268 * t3578;
t3269 = t3539 * t3798 + (t3416 * t3521 + t3671) * t3417 + t3416 * (pkin(1) * t3599 - t3428 * t3575 / 0.2e1) * t3372;
t3524 = -0.2e1 * t3269 * t3576;
t3520 = t3333 * t3431 * t3610;
t3519 = t3334 * t3433 * t3605;
t3518 = t3335 * t3435 * t3600;
t3452 = 1 / pkin(3) ^ 2;
t3517 = t3330 * t3418 * t3452 * t3707;
t3516 = t3331 * t3419 * t3452 * t3701;
t3515 = t3332 * t3420 * t3452 * t3695;
t3240 = t3794 - t3281 / 0.2e1 + (t3644 / 0.2e1 + t3548 / 0.2e1) * t3741;
t3508 = -0.2e1 * t3255 * t3382 - 0.2e1 * t3240;
t3241 = t3793 - t3279 / 0.2e1 + (t3646 / 0.2e1 + t3552 / 0.2e1) * t3751;
t3507 = -0.2e1 * t3256 * t3382 - 0.2e1 * t3241;
t3242 = t3792 - t3280 / 0.2e1 + (t3645 / 0.2e1 + t3550 / 0.2e1) * t3746;
t3506 = -0.2e1 * t3257 * t3382 - 0.2e1 * t3242;
t3505 = t3430 * t3517;
t3504 = t3432 * t3516;
t3503 = t3426 * t3516;
t3502 = t3434 * t3515;
t3501 = t3428 * t3515;
t3500 = t3424 * t3517;
t3487 = (t3343 * t3366 + t3345 * t3705) * t3728;
t3486 = (t3343 * t3705 - t3345 * t3366) * t3728;
t3485 = (t3344 * t3370 + t3346 * t3699) * t3727;
t3484 = (t3344 * t3699 - t3346 * t3370) * t3727;
t3483 = (t3342 * t3374 + t3347 * t3693) * t3726;
t3482 = (t3342 * t3693 - t3347 * t3374) * t3726;
t3480 = (t3337 * t3364 * t3431 * t3663 - t3288 * t3312 * t3365) * t3387;
t3479 = (t3338 * t3368 * t3433 * t3662 - t3289 * t3313 * t3369) * t3390;
t3478 = (t3336 * t3372 * t3435 * t3661 - t3290 * t3314 * t3373) * t3393;
t3477 = 0.2e1 * qJ(2,1) * t3642 + 0.2e1 * t3373 * t3549;
t3476 = 0.2e1 * qJ(2,2) * t3630 + 0.2e1 * t3369 * t3551;
t3475 = 0.2e1 * qJ(2,3) * t3636 + 0.2e1 * t3365 * t3553;
t3247 = -t3473 + 0.2e1 * t3792;
t3245 = -t3474 + 0.2e1 * t3793;
t3243 = -t3472 + 0.2e1 * t3794;
t3239 = (t3449 + t3455) * t3255 - pkin(1) * t3472;
t3238 = (t3448 + t3455) * t3257 - pkin(1) * t3473;
t3237 = (t3447 + t3455) * t3256 - pkin(1) * t3474;
t3236 = 0.2e1 * t3434 * t3655 + (t3243 * t3434 - t3501) * t3417 + (t3428 * t3508 - t3502) * t3416;
t3235 = 0.2e1 * t3432 * t3653 + (t3247 * t3432 - t3503) * t3417 + (t3426 * t3506 - t3504) * t3416;
t3234 = 0.2e1 * t3430 * t3654 + (t3245 * t3430 - t3500) * t3417 + (t3424 * t3507 - t3505) * t3416;
t3233 = t3655 * t3802 + (t3240 * t3802 - t3502) * t3417 + (t3434 * t3508 + t3501) * t3416;
t3232 = t3653 * t3803 + (t3242 * t3803 - t3504) * t3417 + (t3432 * t3506 + t3503) * t3416;
t3231 = t3654 * t3804 + (t3241 * t3804 - t3505) * t3417 + (t3430 * t3507 + t3500) * t3416;
t1 = [(t3255 * t3587 + t3256 * t3591 + t3257 * t3589) * MDP(1) + (-t3274 * t3547 - t3276 * t3546 - t3278 * t3545 + t3322 * t3475 + t3324 * t3476 + t3326 * t3477) * MDP(6) + (t3237 * t3591 + t3238 * t3589 + t3239 * t3587 + (t3276 * t3780 + (-t3276 * t3289 + t3324 * t3662) * t3648) * t3338 + (t3274 * t3781 + (-t3274 * t3288 + t3322 * t3663) * t3647) * t3337 + (t3278 * t3782 + (-t3278 * t3290 + t3326 * t3661) * t3649) * t3336) * MDP(7) + (t3326 * t3556 + t3322 * t3555 + t3324 * t3554 + ((0.4e1 * t3602 - 0.2e1 * t3754) * t3596 + (0.4e1 * t3607 - 0.2e1 * t3756) * t3597 + (0.4e1 * t3612 - 0.2e1 * t3758) * t3598) * t3451) * MDP(8) + (t3326 * t3529 + t3322 * t3528 + t3324 * t3527 + ((0.2e1 * t3602 - t3754) * t3593 + (0.2e1 * t3607 - t3756) * t3594 + (0.2e1 * t3612 - t3758) * t3595) * t3451) * MDP(9) + ((t3345 * t3635 + t3346 * t3629 + t3347 * t3641) * t3451 + (t3322 * t3487 + t3324 * t3485 + t3326 * t3483) * t3452) * MDP(10) + ((t3342 * t3641 + t3343 * t3635 + t3344 * t3629) * t3451 + (t3322 * t3486 + t3324 * t3484 + t3326 * t3482) * t3452) * MDP(11) + (t3394 * t3586 + t3395 * t3584 + t3396 * t3582) * t3783 + ((t3236 * t3729 - t3759) * t3391 + (t3235 * t3731 - t3761) * t3388 + (t3234 * t3733 - t3763) * t3385 + (t3322 * t3526 + t3324 * t3525 + t3326 * t3524 + (-t3536 + (t3278 * t3560 + t3621) * t3740) * t3347 + (-t3537 + (t3276 * t3561 + t3623) * t3745) * t3346 + (-t3538 + (t3274 * t3562 + t3625) * t3750) * t3345) * t3451) * MDP(13) + ((t3233 * t3729 + t3278 * t3643) * t3391 + (t3232 * t3731 + t3276 * t3631) * t3388 + (t3231 * t3733 + t3274 * t3637) * t3385 + (-t3342 * t3536 - t3343 * t3538 - t3344 * t3537 + (t3326 * t3563 + (t3278 * t3658 + t3342 * t3621) * t3293) * t3373 + (t3324 * t3564 + (t3276 * t3659 + t3344 * t3623) * t3292) * t3369 + (t3322 * t3565 + (t3274 * t3660 + t3343 * t3625) * t3291) * t3365) * t3451) * MDP(14) + t3820 * ((t3245 * t3733 - t3337 * t3763) * t3385 + (t3247 * t3731 - t3338 * t3761) * t3388 + (t3243 * t3729 - t3336 * t3759) * t3391); (t3255 * t3588 + t3256 * t3592 + t3257 * t3590) * MDP(1) + (-t3273 * t3547 - t3275 * t3546 - t3277 * t3545 + t3321 * t3475 + t3323 * t3476 + t3325 * t3477) * MDP(6) + (t3237 * t3592 + t3238 * t3590 + t3239 * t3588 + (t3275 * t3780 + (-t3275 * t3289 + t3323 * t3662) * t3648) * t3338 + (t3273 * t3781 + (-t3273 * t3288 + t3321 * t3663) * t3647) * t3337 + (t3277 * t3782 + (-t3277 * t3290 + t3325 * t3661) * t3649) * t3336) * MDP(7) + (t3325 * t3556 + t3321 * t3555 + t3323 * t3554 + ((0.4e1 * t3603 - 0.2e1 * t3753) * t3596 + (0.4e1 * t3608 - 0.2e1 * t3755) * t3597 + (0.4e1 * t3613 - 0.2e1 * t3757) * t3598) * t3451) * MDP(8) + (t3325 * t3529 + t3321 * t3528 + t3323 * t3527 + ((0.2e1 * t3603 - t3753) * t3593 + (0.2e1 * t3608 - t3755) * t3594 + (0.2e1 * t3613 - t3757) * t3595) * t3451) * MDP(9) + ((t3345 * t3634 + t3346 * t3628 + t3347 * t3640) * t3451 + (t3321 * t3487 + t3323 * t3485 + t3325 * t3483) * t3452) * MDP(10) + ((t3342 * t3640 + t3343 * t3634 + t3344 * t3628) * t3451 + (t3321 * t3486 + t3323 * t3484 + t3325 * t3482) * t3452) * MDP(11) + (t3397 * t3586 + t3398 * t3584 + t3399 * t3582) * t3783 + ((t3236 * t3730 - t3760) * t3391 + (t3235 * t3732 - t3762) * t3388 + (t3234 * t3734 - t3764) * t3385 + (t3321 * t3526 + t3323 * t3525 + t3325 * t3524 + (-t3533 + (t3277 * t3560 + t3620) * t3740) * t3347 + (-t3534 + (t3275 * t3561 + t3622) * t3745) * t3346 + (-t3535 + (t3273 * t3562 + t3624) * t3750) * t3345) * t3451) * MDP(13) + ((t3233 * t3730 + t3277 * t3643) * t3391 + (t3232 * t3732 + t3275 * t3631) * t3388 + (t3231 * t3734 + t3273 * t3637) * t3385 + (-t3342 * t3533 - t3343 * t3535 - t3344 * t3534 + (t3325 * t3563 + (t3277 * t3658 + t3342 * t3620) * t3293) * t3373 + (t3323 * t3564 + (t3275 * t3659 + t3344 * t3622) * t3292) * t3369 + (t3321 * t3565 + (t3273 * t3660 + t3343 * t3624) * t3291) * t3365) * t3451) * MDP(14) + t3820 * (t3385 * (t3245 * t3734 - t3337 * t3764) + t3388 * (t3247 * t3732 - t3338 * t3762) + t3391 * (t3243 * t3730 - t3336 * t3760)); (t3626 + t3632 + t3638) * MDP(1) + (qJ(2,1) * t3568 + qJ(2,2) * t3566 + qJ(2,3) * t3567 + t3478 + t3479 + t3480) * MDP(6) + ((t3435 * t3239 + t3314 * t3244) * t3391 + (t3433 * t3238 + t3313 * t3248) * t3388 + (t3431 * t3237 + t3312 * t3246) * t3385 + qJ(2,3) * t3480 + qJ(2,2) * t3479 + qJ(2,1) * t3478) * MDP(7) + (t3341 * t3638 + t3339 * t3632 + t3340 * t3626 + 0.4e1 * (t3303 * t3520 + t3304 * t3519 + t3305 * t3518) * t3451) * MDP(8) + 0.2e1 * (t3305 * t3568 + t3303 * t3567 + t3304 * t3566 + (t3309 * t3520 + t3310 * t3519 + t3311 * t3518) * t3451) * MDP(9) + ((t3236 * t3435 - t3342 * t3779) * t3391 + (t3235 * t3433 - t3344 * t3773) * t3388 + (t3234 * t3431 - t3343 * t3776) * t3385 + ((-t3269 * t3572 + t3347 * t3540) * t3817 + (-t3268 * t3573 + t3346 * t3542) * t3818 + (-t3267 * t3574 + t3345 * t3544) * t3819) * t3451) * MDP(13) + ((t3233 * t3435 + t3314 * t3778) * t3391 + (t3232 * t3433 + t3313 * t3772) * t3388 + (t3231 * t3431 + t3312 * t3775) * t3385 + ((t3272 * t3572 + t3342 * t3540) * t3817 + (t3271 * t3573 + t3344 * t3542) * t3818 + (t3270 * t3574 + t3343 * t3544) * t3819) * t3451) * MDP(14) + (((t3342 * t3373 + t3347 * t3695) * t3581 + (t3344 * t3369 + t3346 * t3701) * t3583 + (t3343 * t3365 + t3345 * t3707) * t3585) * MDP(10) + ((t3342 * t3695 - t3347 * t3373) * t3581 + (t3344 * t3701 - t3346 * t3369) * t3583 + (t3343 * t3707 - t3345 * t3365) * t3585) * MDP(11)) * t3452 + t3820 * (t3385 * (t3245 * t3431 - t3776) + t3388 * (t3247 * t3433 - t3773) + t3391 * (t3243 * t3435 - t3779));];
taucX  = t1;