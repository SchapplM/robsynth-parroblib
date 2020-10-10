% Calculate Gravitation load for parallel robot
% P6RRPRRR14V3G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-15 09:53
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(1,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V3G1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-15 09:48:59
% EndTime: 2019-04-15 09:49:04
% DurationCPUTime: 5.43s
% Computational Cost: add. (2141->293), mult. (4125->525), div. (126->12), fcn. (2783->42), ass. (0->233)
t3584 = xP(5);
t3535 = sin(t3584);
t3538 = cos(t3584);
t3595 = koppelP(6,3);
t3583 = xP(6);
t3534 = sin(t3583);
t3537 = cos(t3583);
t3601 = koppelP(6,2);
t3607 = koppelP(6,1);
t3710 = -t3534 * t3601 + t3537 * t3607;
t3453 = -t3535 * t3710 + t3538 * t3595;
t3484 = t3534 * t3607 + t3537 * t3601;
t3585 = xP(4);
t3536 = sin(t3585);
t3539 = cos(t3585);
t3417 = t3453 * t3539 + t3484 * t3536;
t3589 = 0.1e1 / qJ(3,6);
t3552 = legFrame(6,3);
t3522 = sin(t3552);
t3528 = cos(t3552);
t3490 = -g(1) * t3522 + g(2) * t3528;
t3496 = g(1) * t3528 + g(2) * t3522;
t3508 = m(1) * rSges(1,2) - m(2) * rSges(2,3) - m(3) * rSges(3,2);
t3558 = sin(qJ(2,6));
t3559 = sin(qJ(1,6));
t3565 = cos(qJ(1,6));
t3582 = m(2) * rSges(2,2);
t3509 = (rSges(3,3) + qJ(3,6)) * m(3) - t3582;
t3521 = m(2) * rSges(2,1) + m(3) * rSges(3,1);
t3564 = cos(qJ(2,6));
t3709 = m(1) * rSges(1,1);
t3618 = t3509 * t3558 + t3521 * t3564 + t3709;
t3701 = ((t3508 * t3565 + t3559 * t3618) * t3496 + (t3508 * t3559 - t3565 * t3618) * t3490) / t3558;
t3637 = t3589 * t3701;
t3722 = t3417 * t3637;
t3596 = koppelP(5,3);
t3602 = koppelP(5,2);
t3608 = koppelP(5,1);
t3711 = -t3534 * t3602 + t3537 * t3608;
t3455 = -t3535 * t3711 + t3538 * t3596;
t3485 = t3534 * t3608 + t3537 * t3602;
t3418 = t3455 * t3539 + t3485 * t3536;
t3590 = 0.1e1 / qJ(3,5);
t3553 = legFrame(5,3);
t3523 = sin(t3553);
t3529 = cos(t3553);
t3491 = -g(1) * t3523 + g(2) * t3529;
t3497 = g(1) * t3529 + g(2) * t3523;
t3560 = sin(qJ(2,5));
t3561 = sin(qJ(1,5));
t3567 = cos(qJ(1,5));
t3510 = (rSges(3,3) + qJ(3,5)) * m(3) - t3582;
t3566 = cos(qJ(2,5));
t3617 = t3510 * t3560 + t3521 * t3566 + t3709;
t3700 = ((t3508 * t3567 + t3561 * t3617) * t3497 + (t3508 * t3561 - t3567 * t3617) * t3491) / t3560;
t3636 = t3590 * t3700;
t3721 = t3418 * t3636;
t3597 = koppelP(4,3);
t3603 = koppelP(4,2);
t3609 = koppelP(4,1);
t3712 = -t3534 * t3603 + t3537 * t3609;
t3457 = -t3535 * t3712 + t3538 * t3597;
t3486 = t3534 * t3609 + t3537 * t3603;
t3419 = t3457 * t3539 + t3486 * t3536;
t3591 = 0.1e1 / qJ(3,4);
t3554 = legFrame(4,3);
t3524 = sin(t3554);
t3530 = cos(t3554);
t3492 = -g(1) * t3524 + g(2) * t3530;
t3498 = g(1) * t3530 + g(2) * t3524;
t3562 = sin(qJ(2,4));
t3563 = sin(qJ(1,4));
t3569 = cos(qJ(1,4));
t3511 = (rSges(3,3) + qJ(3,4)) * m(3) - t3582;
t3568 = cos(qJ(2,4));
t3616 = t3511 * t3562 + t3521 * t3568 + t3709;
t3699 = ((t3508 * t3569 + t3563 * t3616) * t3498 + (t3508 * t3563 - t3569 * t3616) * t3492) / t3562;
t3635 = t3591 * t3699;
t3720 = t3419 * t3635;
t3598 = koppelP(3,3);
t3604 = koppelP(3,2);
t3610 = koppelP(3,1);
t3713 = -t3534 * t3604 + t3537 * t3610;
t3459 = -t3535 * t3713 + t3538 * t3598;
t3487 = t3534 * t3610 + t3537 * t3604;
t3420 = t3459 * t3539 + t3487 * t3536;
t3592 = 0.1e1 / qJ(3,3);
t3555 = legFrame(3,3);
t3525 = sin(t3555);
t3531 = cos(t3555);
t3493 = -g(1) * t3525 + g(2) * t3531;
t3499 = g(1) * t3531 + g(2) * t3525;
t3570 = sin(qJ(2,3));
t3571 = sin(qJ(1,3));
t3577 = cos(qJ(1,3));
t3515 = (rSges(3,3) + qJ(3,3)) * m(3) - t3582;
t3576 = cos(qJ(2,3));
t3615 = t3515 * t3570 + t3521 * t3576 + t3709;
t3698 = ((t3508 * t3577 + t3571 * t3615) * t3499 + (t3508 * t3571 - t3577 * t3615) * t3493) / t3570;
t3634 = t3592 * t3698;
t3719 = t3420 * t3634;
t3599 = koppelP(2,3);
t3605 = koppelP(2,2);
t3611 = koppelP(2,1);
t3714 = -t3534 * t3605 + t3537 * t3611;
t3461 = -t3535 * t3714 + t3538 * t3599;
t3488 = t3534 * t3611 + t3537 * t3605;
t3434 = t3461 * t3539 + t3488 * t3536;
t3593 = 0.1e1 / qJ(3,2);
t3556 = legFrame(2,3);
t3526 = sin(t3556);
t3532 = cos(t3556);
t3494 = -g(1) * t3526 + g(2) * t3532;
t3500 = g(1) * t3532 + g(2) * t3526;
t3572 = sin(qJ(2,2));
t3573 = sin(qJ(1,2));
t3579 = cos(qJ(1,2));
t3516 = (rSges(3,3) + qJ(3,2)) * m(3) - t3582;
t3578 = cos(qJ(2,2));
t3614 = t3516 * t3572 + t3521 * t3578 + t3709;
t3697 = ((t3508 * t3579 + t3573 * t3614) * t3500 + (t3508 * t3573 - t3579 * t3614) * t3494) / t3572;
t3633 = t3593 * t3697;
t3718 = t3434 * t3633;
t3600 = koppelP(1,3);
t3606 = koppelP(1,2);
t3612 = koppelP(1,1);
t3715 = -t3534 * t3606 + t3537 * t3612;
t3463 = -t3535 * t3715 + t3538 * t3600;
t3489 = t3534 * t3612 + t3537 * t3606;
t3438 = t3463 * t3539 + t3489 * t3536;
t3594 = 0.1e1 / qJ(3,1);
t3557 = legFrame(1,3);
t3527 = sin(t3557);
t3533 = cos(t3557);
t3495 = -g(1) * t3527 + g(2) * t3533;
t3501 = g(1) * t3533 + g(2) * t3527;
t3574 = sin(qJ(2,1));
t3575 = sin(qJ(1,1));
t3581 = cos(qJ(1,1));
t3517 = (rSges(3,3) + qJ(3,1)) * m(3) - t3582;
t3580 = cos(qJ(2,1));
t3613 = t3517 * t3574 + t3521 * t3580 + t3709;
t3696 = ((t3508 * t3581 + t3575 * t3613) * t3501 + (t3508 * t3575 - t3581 * t3613) * t3495) / t3574;
t3632 = t3594 * t3696;
t3717 = t3438 * t3632;
t3421 = t3453 * t3536 - t3484 * t3539;
t3424 = t3455 * t3536 - t3485 * t3539;
t3427 = t3457 * t3536 - t3486 * t3539;
t3430 = t3459 * t3536 - t3487 * t3539;
t3433 = t3461 * t3536 - t3488 * t3539;
t3437 = t3463 * t3536 - t3489 * t3539;
t3586 = rSges(4,3);
t3587 = rSges(4,2);
t3588 = rSges(4,1);
t3619 = t3534 * t3587 - t3537 * t3588;
t3716 = -t3535 * t3586 + t3538 * t3619;
t3625 = t3490 * t3559 + t3496 * t3565;
t3441 = -g(3) * t3564 + t3558 * t3625;
t3708 = m(3) * t3441;
t3624 = t3491 * t3561 + t3497 * t3567;
t3442 = -g(3) * t3566 + t3560 * t3624;
t3707 = m(3) * t3442;
t3623 = t3492 * t3563 + t3498 * t3569;
t3443 = -g(3) * t3568 + t3562 * t3623;
t3706 = m(3) * t3443;
t3622 = t3493 * t3571 + t3499 * t3577;
t3444 = -g(3) * t3576 + t3570 * t3622;
t3705 = m(3) * t3444;
t3621 = t3494 * t3573 + t3500 * t3579;
t3445 = -g(3) * t3578 + t3572 * t3621;
t3704 = m(3) * t3445;
t3620 = t3495 * t3575 + t3501 * t3581;
t3446 = -g(3) * t3580 + t3574 * t3620;
t3703 = m(3) * t3446;
t3702 = g(3) * t3521;
t3411 = (-t3509 * t3625 - t3702) * t3564 + t3558 * (-g(3) * t3509 + t3521 * t3625);
t3695 = t3411 * t3589;
t3412 = (-t3510 * t3624 - t3702) * t3566 + t3560 * (-g(3) * t3510 + t3521 * t3624);
t3694 = t3412 * t3590;
t3413 = (-t3511 * t3623 - t3702) * t3568 + t3562 * (-g(3) * t3511 + t3521 * t3623);
t3693 = t3413 * t3591;
t3414 = (-t3515 * t3622 - t3702) * t3576 + t3570 * (-g(3) * t3515 + t3521 * t3622);
t3692 = t3414 * t3592;
t3415 = (-t3516 * t3621 - t3702) * t3578 + t3572 * (-g(3) * t3516 + t3521 * t3621);
t3691 = t3415 * t3593;
t3416 = (-t3517 * t3620 - t3702) * t3580 + t3574 * (-g(3) * t3517 + t3521 * t3620);
t3690 = t3416 * t3594;
t3465 = -t3522 * t3559 + t3528 * t3565;
t3689 = t3465 * t3558;
t3688 = t3465 * t3564;
t3466 = t3522 * t3565 + t3528 * t3559;
t3687 = t3466 * t3558;
t3686 = t3466 * t3564;
t3467 = -t3523 * t3561 + t3529 * t3567;
t3685 = t3467 * t3560;
t3684 = t3467 * t3566;
t3468 = t3523 * t3567 + t3529 * t3561;
t3683 = t3468 * t3560;
t3682 = t3468 * t3566;
t3469 = -t3524 * t3563 + t3530 * t3569;
t3681 = t3469 * t3562;
t3680 = t3469 * t3568;
t3470 = t3524 * t3569 + t3530 * t3563;
t3679 = t3470 * t3562;
t3678 = t3470 * t3568;
t3471 = -t3525 * t3571 + t3531 * t3577;
t3677 = t3471 * t3570;
t3676 = t3471 * t3576;
t3472 = t3525 * t3577 + t3531 * t3571;
t3675 = t3472 * t3570;
t3674 = t3472 * t3576;
t3473 = -t3526 * t3573 + t3532 * t3579;
t3673 = t3473 * t3572;
t3672 = t3473 * t3578;
t3474 = t3526 * t3579 + t3532 * t3573;
t3671 = t3474 * t3572;
t3670 = t3474 * t3578;
t3475 = -t3527 * t3575 + t3533 * t3581;
t3669 = t3475 * t3574;
t3668 = t3475 * t3580;
t3476 = t3527 * t3581 + t3533 * t3575;
t3667 = t3476 * t3574;
t3666 = t3476 * t3580;
t3646 = t3535 * t3587;
t3645 = t3535 * t3588;
t3638 = t3538 * t3586;
t3452 = t3535 * t3600 + t3538 * t3715;
t3451 = t3535 * t3599 + t3538 * t3714;
t3450 = t3535 * t3598 + t3538 * t3713;
t3449 = t3535 * t3597 + t3538 * t3712;
t3448 = t3535 * t3596 + t3538 * t3711;
t3447 = t3535 * t3595 + t3538 * t3710;
t1 = [-m(4) * g(1) + (t3416 * t3668 - t3476 * t3696) * t3594 + (t3415 * t3672 - t3474 * t3697) * t3593 + (t3414 * t3676 - t3472 * t3698) * t3592 + (t3413 * t3680 - t3470 * t3699) * t3591 + (t3412 * t3684 - t3468 * t3700) * t3590 + (t3411 * t3688 - t3466 * t3701) * t3589 + (-t3441 * t3689 - t3442 * t3685 - t3443 * t3681 - t3444 * t3677 - t3445 * t3673 - t3446 * t3669) * m(3); -m(4) * g(2) + (t3416 * t3666 + t3475 * t3696) * t3594 + (t3415 * t3670 + t3473 * t3697) * t3593 + (t3414 * t3674 + t3471 * t3698) * t3592 + (t3413 * t3678 + t3469 * t3699) * t3591 + (t3412 * t3682 + t3467 * t3700) * t3590 + (t3411 * t3686 + t3465 * t3701) * t3589 + (-t3441 * t3687 - t3442 * t3683 - t3443 * t3679 - t3444 * t3675 - t3445 * t3671 - t3446 * t3667) * m(3); t3558 * t3695 + t3560 * t3694 + t3562 * t3693 + t3570 * t3692 + t3572 * t3691 + t3574 * t3690 - m(4) * g(3) + (t3441 * t3564 + t3442 * t3566 + t3443 * t3568 + t3444 * t3576 + t3445 * t3578 + t3446 * t3580) * m(3); -t3475 * t3717 + (-t3437 * t3574 - t3438 * t3666) * t3690 - (t3437 * t3580 - t3438 * t3667) * t3703 - t3473 * t3718 + (-t3433 * t3572 - t3434 * t3670) * t3691 - (t3433 * t3578 - t3434 * t3671) * t3704 - t3471 * t3719 + (-t3420 * t3674 - t3430 * t3570) * t3692 - (-t3420 * t3675 + t3430 * t3576) * t3705 - t3469 * t3720 + (-t3419 * t3678 - t3427 * t3562) * t3693 - (-t3419 * t3679 + t3427 * t3568) * t3706 - t3467 * t3721 + (-t3418 * t3682 - t3424 * t3560) * t3694 - (-t3418 * t3683 + t3424 * t3566) * t3707 - t3465 * t3722 + (-t3417 * t3686 - t3421 * t3558) * t3695 - (-t3417 * t3687 + t3421 * t3564) * t3708 + m(4) * (((-g(2) * t3645 - g(3) * t3587) * t3537 + g(2) * t3638) * t3539 + ((g(2) * t3587 - g(3) * t3645) * t3537 + g(3) * t3638) * t3536 + ((g(2) * t3646 - g(3) * t3588) * t3539 + (g(2) * t3588 + g(3) * t3646) * t3536) * t3534); -t3476 * t3717 + (t3438 * t3668 - t3452 * t3574) * t3690 - (t3438 * t3669 + t3452 * t3580) * t3703 - t3474 * t3718 + (t3434 * t3672 - t3451 * t3572) * t3691 - (t3434 * t3673 + t3451 * t3578) * t3704 - t3472 * t3719 + (t3420 * t3676 - t3450 * t3570) * t3692 - (t3420 * t3677 + t3450 * t3576) * t3705 - t3470 * t3720 + (t3419 * t3680 - t3449 * t3562) * t3693 - (t3419 * t3681 + t3449 * t3568) * t3706 - t3468 * t3721 + (t3418 * t3684 - t3448 * t3560) * t3694 - (t3418 * t3685 + t3448 * t3566) * t3707 - t3466 * t3722 + (t3417 * t3688 - t3447 * t3558) * t3695 - (t3417 * t3689 + t3447 * t3564) * t3708 - m(4) * (t3716 * g(3) + ((t3535 * t3619 + t3638) * t3539 + (t3534 * t3588 + t3537 * t3587) * t3536) * g(1)); (-t3437 * t3476 + t3452 * t3475) * t3632 + (-t3433 * t3474 + t3451 * t3473) * t3633 + (-t3430 * t3472 + t3450 * t3471) * t3634 + (-t3427 * t3470 + t3449 * t3469) * t3635 + (-t3424 * t3468 + t3448 * t3467) * t3636 + (-t3421 * t3466 + t3447 * t3465) * t3637 + m(4) * (t3716 * g(2) + (-t3536 * t3638 + (t3536 * t3645 + t3539 * t3587) * t3537 + (-t3536 * t3646 + t3539 * t3588) * t3534) * g(1)) + (-t3558 * t3708 + t3564 * t3695) * (t3421 * t3465 + t3447 * t3466) + (-t3560 * t3707 + t3566 * t3694) * (t3424 * t3467 + t3448 * t3468) + (-t3562 * t3706 + t3568 * t3693) * (t3427 * t3469 + t3449 * t3470) + (-t3570 * t3705 + t3576 * t3692) * (t3430 * t3471 + t3450 * t3472) + (-t3572 * t3704 + t3578 * t3691) * (t3433 * t3473 + t3451 * t3474) + (-t3574 * t3703 + t3580 * t3690) * (t3437 * t3475 + t3452 * t3476);];
taugX  = t1;
