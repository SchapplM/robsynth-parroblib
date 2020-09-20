% Calculate inertia matrix for parallel robot
% P3RRRRR2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:10
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:08:37
% EndTime: 2020-03-09 21:08:38
% DurationCPUTime: 1.30s
% Computational Cost: add. (3168->222), mult. (6531->405), div. (759->14), fcn. (2835->45), ass. (0->216)
t594 = cos(qJ(3,1));
t718 = t594 ^ 2;
t591 = cos(qJ(3,2));
t717 = t591 ^ 2;
t588 = cos(qJ(3,3));
t716 = t588 ^ 2;
t710 = m(3) / 0.2e1;
t715 = Icges(3,2) / 0.2e1;
t709 = m(3) * rSges(3,2);
t553 = -rSges(3,1) * t709 + Icges(3,4);
t603 = 0.2e1 * qJ(3,3);
t606 = rSges(3,2) ^ 2;
t607 = rSges(3,1) ^ 2;
t707 = (-t606 + t607) * t710 - Icges(3,1) / 0.2e1 + t715;
t714 = cos(t603) * t707 + t553 * sin(t603);
t604 = 0.2e1 * qJ(3,2);
t713 = cos(t604) * t707 + t553 * sin(t604);
t605 = 0.2e1 * qJ(3,1);
t712 = cos(t605) * t707 + t553 * sin(t605);
t711 = -2 * pkin(1);
t708 = m(3) * rSges(3,3);
t581 = sin(qJ(1,3));
t706 = pkin(1) * t581;
t584 = sin(qJ(1,2));
t705 = pkin(1) * t584;
t587 = sin(qJ(1,1));
t704 = pkin(1) * t587;
t576 = legFrame(3,2);
t557 = sin(t576);
t560 = cos(t576);
t579 = sin(qJ(3,3));
t681 = t560 * t579;
t580 = sin(qJ(2,3));
t589 = cos(qJ(2,3));
t590 = cos(qJ(1,3));
t538 = t590 * t580 + t581 * t589;
t694 = t538 * t588;
t517 = -t557 * t694 + t681;
t568 = 0.1e1 / t588;
t703 = t517 * t568;
t687 = t557 * t579;
t518 = t560 * t694 + t687;
t702 = t518 * t568;
t577 = legFrame(2,2);
t558 = sin(t577);
t561 = cos(t577);
t582 = sin(qJ(3,2));
t679 = t561 * t582;
t583 = sin(qJ(2,2));
t592 = cos(qJ(2,2));
t593 = cos(qJ(1,2));
t539 = t593 * t583 + t584 * t592;
t693 = t539 * t591;
t519 = -t558 * t693 + t679;
t571 = 0.1e1 / t591;
t701 = t519 * t571;
t685 = t558 * t582;
t520 = t561 * t693 + t685;
t700 = t520 * t571;
t578 = legFrame(1,2);
t559 = sin(t578);
t562 = cos(t578);
t585 = sin(qJ(3,1));
t677 = t562 * t585;
t586 = sin(qJ(2,1));
t595 = cos(qJ(2,1));
t596 = cos(qJ(1,1));
t540 = t596 * t586 + t587 * t595;
t692 = t540 * t594;
t521 = -t559 * t692 + t677;
t574 = 0.1e1 / t594;
t699 = t521 * t574;
t683 = t559 * t585;
t522 = t562 * t692 + t683;
t698 = t522 * t574;
t656 = qJ(2,3) + qJ(3,3);
t657 = qJ(2,3) - qJ(3,3);
t697 = (t590 * t711 + (-cos(qJ(1,3) + t657) - cos(qJ(1,3) + t656)) * pkin(2)) / (sin(t656) + sin(t657));
t658 = qJ(2,2) + qJ(3,2);
t659 = qJ(2,2) - qJ(3,2);
t696 = (t593 * t711 + (-cos(qJ(1,2) + t659) - cos(qJ(1,2) + t658)) * pkin(2)) / (sin(t658) + sin(t659));
t660 = qJ(2,1) + qJ(3,1);
t661 = qJ(2,1) - qJ(3,1);
t695 = (t596 * t711 + (-cos(qJ(1,1) + t661) - cos(qJ(1,1) + t660)) * pkin(2)) / (sin(t660) + sin(t661));
t564 = 0.1e1 / t580;
t691 = cos(qJ(1,3) + qJ(2,3)) * t564;
t565 = 0.1e1 / t583;
t690 = cos(qJ(1,2) + qJ(2,2)) * t565;
t566 = 0.1e1 / t586;
t689 = cos(qJ(1,1) + qJ(2,1)) * t566;
t688 = t557 * t568;
t686 = t558 * t571;
t684 = t559 * t574;
t682 = t560 * t568;
t680 = t561 * t571;
t678 = t562 * t574;
t676 = t564 * t568;
t569 = 0.1e1 / t716;
t675 = t564 * t569;
t610 = 1 / pkin(1);
t674 = t564 * t610;
t673 = t565 * t571;
t572 = 0.1e1 / t717;
t672 = t565 * t572;
t671 = t565 * t610;
t670 = t566 * t574;
t575 = 0.1e1 / t718;
t669 = t566 * t575;
t668 = t566 * t610;
t608 = 0.1e1 / pkin(2);
t667 = t568 * t608;
t666 = t569 * t608;
t665 = t571 * t608;
t664 = t572 * t608;
t663 = t574 * t608;
t662 = t575 * t608;
t552 = rSges(3,1) * t708 - Icges(3,5);
t655 = t606 + t607;
t654 = m(3) * rSges(3,1) * pkin(1);
t653 = t538 * t716 * pkin(2);
t652 = t539 * t717 * pkin(2);
t651 = t540 * t718 * pkin(2);
t650 = t589 * t579 * pkin(1);
t649 = t592 * t582 * pkin(1);
t648 = t595 * t585 * pkin(1);
t616 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + t715 + Icges(3,1) / 0.2e1;
t617 = 0.2e1 * rSges(3,3) ^ 2 + t655;
t615 = t617 * t710 + t616;
t511 = t615 + t714;
t550 = m(2) * rSges(2,2) - t708;
t600 = m(2) * rSges(2,1);
t614 = (-(-t600 + (-rSges(3,1) * t588 + rSges(3,2) * t579) * m(3)) * t589 - t550 * t580) * pkin(1);
t502 = t614 + t511;
t647 = t502 * t666;
t512 = t615 + t713;
t613 = (-(-t600 + (-rSges(3,1) * t591 + rSges(3,2) * t582) * m(3)) * t592 - t550 * t583) * pkin(1);
t503 = t613 + t512;
t646 = t503 * t664;
t513 = t615 + t712;
t612 = (-(-t600 + (-rSges(3,1) * t594 + rSges(3,2) * t585) * m(3)) * t595 - t550 * t586) * pkin(1);
t504 = t612 + t513;
t645 = t504 * t662;
t505 = t557 * t653 + (-pkin(2) * t681 + t557 * t706) * t588 - t560 * t650;
t644 = t505 * t675;
t506 = -t560 * t653 + (-pkin(2) * t687 - t560 * t706) * t588 - t557 * t650;
t643 = t506 * t675;
t507 = t558 * t652 + (-pkin(2) * t679 + t558 * t705) * t591 - t561 * t649;
t642 = t507 * t672;
t508 = -t561 * t652 + (-pkin(2) * t685 - t561 * t705) * t591 - t558 * t649;
t641 = t508 * t672;
t509 = t559 * t651 + (-pkin(2) * t677 + t559 * t704) * t594 - t562 * t648;
t640 = t509 * t669;
t510 = -t562 * t651 + (-pkin(2) * t683 - t562 * t704) * t594 - t559 * t648;
t639 = t510 * t669;
t638 = t511 * t666;
t637 = t512 * t664;
t636 = t513 * t662;
t635 = t517 * t676;
t634 = t518 * t676;
t633 = t519 * t673;
t632 = t520 * t673;
t631 = t521 * t670;
t630 = t522 * t670;
t629 = t608 * t697;
t628 = t608 * t696;
t627 = t608 * t695;
t551 = -rSges(3,2) * t708 + Icges(3,6);
t526 = t551 * t588 - t552 * t579;
t626 = t526 * t666;
t527 = t551 * t591 - t552 * t582;
t625 = t527 * t664;
t528 = t551 * t594 - t552 * t585;
t624 = t528 * t662;
t623 = t557 * t667;
t622 = t558 * t665;
t621 = t559 * t663;
t620 = t560 * t667;
t619 = t561 * t665;
t618 = t562 * t663;
t609 = pkin(1) ^ 2;
t611 = Icges(1,3) + ((2 * t609) + t617) * t710 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + m(2) * t609 + t616;
t549 = t655 * m(3) + Icges(3,3);
t516 = (Icges(3,6) + (-pkin(1) * t586 - rSges(3,3)) * t709) * t594 - t585 * (t586 * t654 + t552);
t515 = (Icges(3,6) + (-pkin(1) * t583 - rSges(3,3)) * t709) * t591 - t582 * (t583 * t654 + t552);
t514 = (Icges(3,6) + (-pkin(1) * t580 - rSges(3,3)) * t709) * t588 - t579 * (t580 * t654 + t552);
t501 = t611 + 0.2e1 * t612 + t712;
t500 = t611 + 0.2e1 * t613 + t713;
t499 = t611 + 0.2e1 * t614 + t714;
t498 = (t516 * t689 + t528 * t627) * t610;
t497 = (t515 * t690 + t527 * t628) * t610;
t496 = (t514 * t691 + t526 * t629) * t610;
t495 = t549 * t621 + (t510 * t624 + t516 * t698) * t668;
t494 = t549 * t618 + (t509 * t624 + t516 * t699) * t668;
t493 = t549 * t622 + (t508 * t625 + t515 * t700) * t671;
t492 = t549 * t619 + (t507 * t625 + t515 * t701) * t671;
t491 = t549 * t623 + (t506 * t626 + t514 * t702) * t674;
t490 = t549 * t620 + (t505 * t626 + t514 * t703) * t674;
t489 = (t504 * t689 + t513 * t627) * t610;
t488 = (t503 * t690 + t512 * t628) * t610;
t487 = (t502 * t691 + t511 * t629) * t610;
t486 = (t501 * t689 + t504 * t627) * t610;
t485 = (t500 * t690 + t503 * t628) * t610;
t484 = (t499 * t691 + t502 * t629) * t610;
t483 = t528 * t621 + (t504 * t698 + t510 * t636) * t668;
t482 = t528 * t618 + (t504 * t699 + t509 * t636) * t668;
t481 = t527 * t622 + (t503 * t700 + t508 * t637) * t671;
t480 = t527 * t619 + (t503 * t701 + t507 * t637) * t671;
t479 = t526 * t623 + (t502 * t702 + t506 * t638) * t674;
t478 = t526 * t620 + (t502 * t703 + t505 * t638) * t674;
t477 = t516 * t621 + (t501 * t698 + t510 * t645) * t668;
t476 = t516 * t618 + (t501 * t699 + t509 * t645) * t668;
t475 = t515 * t622 + (t500 * t700 + t508 * t646) * t671;
t474 = t515 * t619 + (t500 * t701 + t507 * t646) * t671;
t473 = t514 * t623 + (t499 * t702 + t506 * t647) * t674;
t472 = t514 * t620 + (t499 * t703 + t505 * t647) * t674;
t1 = [m(4) + (t473 * t634 + t475 * t632 + t477 * t630) * t610 + (t491 * t688 + t493 * t686 + t495 * t684 + (t479 * t643 + t481 * t641 + t483 * t639) * t610) * t608, (t473 * t635 + t475 * t633 + t477 * t631) * t610 + (t491 * t682 + t493 * t680 + t495 * t678 + (t479 * t644 + t481 * t642 + t483 * t640) * t610) * t608, (t473 * t691 + t475 * t690 + t477 * t689 + (t479 * t697 + t481 * t696 + t483 * t695) * t608) * t610; (t472 * t634 + t474 * t632 + t476 * t630) * t610 + (t490 * t688 + t492 * t686 + t494 * t684 + (t478 * t643 + t480 * t641 + t482 * t639) * t610) * t608, m(4) + (t472 * t635 + t474 * t633 + t476 * t631) * t610 + (t490 * t682 + t492 * t680 + t494 * t678 + (t478 * t644 + t480 * t642 + t482 * t640) * t610) * t608, (t472 * t691 + t474 * t690 + t476 * t689 + (t478 * t697 + t480 * t696 + t482 * t695) * t608) * t610; (t484 * t634 + t485 * t632 + t486 * t630) * t610 + (t496 * t688 + t497 * t686 + t498 * t684 + (t487 * t643 + t488 * t641 + t489 * t639) * t610) * t608, (t484 * t635 + t485 * t633 + t486 * t631) * t610 + (t496 * t682 + t497 * t680 + t498 * t678 + (t487 * t644 + t488 * t642 + t489 * t640) * t610) * t608, m(4) + (t484 * t691 + t485 * t690 + t486 * t689 + (t487 * t697 + t488 * t696 + t489 * t695) * t608) * t610;];
MX  = t1;
