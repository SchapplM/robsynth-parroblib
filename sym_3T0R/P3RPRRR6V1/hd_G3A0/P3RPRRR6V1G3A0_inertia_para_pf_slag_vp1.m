% Calculate inertia matrix for parallel robot
% P3RPRRR6V1G3A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:18
% EndTime: 2020-08-06 18:41:21
% DurationCPUTime: 1.99s
% Computational Cost: add. (4206->263), mult. (6048->448), div. (531->10), fcn. (3690->56), ass. (0->185)
t621 = sin(pkin(7));
t642 = -pkin(6) - pkin(5);
t583 = t642 * t621 - pkin(1);
t633 = cos(qJ(1,3));
t567 = t583 * t633;
t622 = cos(pkin(7));
t627 = sin(qJ(1,3));
t694 = t627 * t642;
t695 = t627 * t621;
t735 = t567 + pkin(2) * t695 - (pkin(2) * t633 - t694) * t622;
t635 = cos(qJ(1,2));
t568 = t583 * t635;
t629 = sin(qJ(1,2));
t690 = t629 * t642;
t691 = t629 * t621;
t734 = t568 + pkin(2) * t691 - (pkin(2) * t635 - t690) * t622;
t637 = cos(qJ(1,1));
t569 = t583 * t637;
t631 = sin(qJ(1,1));
t686 = t631 * t642;
t687 = t631 * t621;
t733 = t569 + pkin(2) * t687 - (pkin(2) * t637 - t686) * t622;
t684 = 2 * m(3);
t732 = m(3) / 0.2e1;
t731 = Icges(3,2) / 0.2e1;
t630 = sin(qJ(3,1));
t636 = cos(qJ(3,1));
t578 = t636 * rSges(3,1) - t630 * rSges(3,2);
t628 = sin(qJ(3,2));
t634 = cos(qJ(3,2));
t577 = t634 * rSges(3,1) - t628 * rSges(3,2);
t626 = sin(qJ(3,3));
t632 = cos(qJ(3,3));
t576 = t632 * rSges(3,1) - t626 * rSges(3,2);
t730 = 0.2e1 * pkin(1);
t729 = 0.2e1 * pkin(2);
t728 = 0.2e1 * t642;
t727 = m(3) * rSges(3,2);
t648 = rSges(3,2) ^ 2;
t649 = rSges(3,1) ^ 2;
t726 = (-t648 + t649) * t732 + t731 - Icges(3,1) / 0.2e1;
t725 = m(3) * t576;
t724 = m(3) * t577;
t723 = m(3) * t578;
t722 = pkin(1) * t621;
t599 = t622 * pkin(1);
t592 = t632 * pkin(3) + pkin(2);
t593 = t634 * pkin(3) + pkin(2);
t594 = t636 * pkin(3) + pkin(2);
t650 = 0.1e1 / pkin(3);
t700 = t592 * t633;
t715 = ((-t694 + t700) * t622 - t567 - t592 * t695) * t650;
t699 = t593 * t635;
t714 = ((-t690 + t699) * t622 - t568 - t593 * t691) * t650;
t698 = t594 * t637;
t713 = ((-t686 + t698) * t622 - t569 - t594 * t687) * t650;
t609 = qJ(1,3) + pkin(7);
t596 = cos(t609);
t645 = 0.2e1 * qJ(3,3);
t612 = sin(t645);
t678 = pkin(7) + qJ(3,3);
t681 = -pkin(7) + qJ(3,3);
t712 = (t596 * t728 + sin(t609) * t729 + t627 * t730 + (sin(qJ(1,3) - t681) + sin(qJ(1,3) + t678)) * pkin(3)) / (t626 * t729 + pkin(3) * t612 + (sin(t678) + sin(t681)) * pkin(1));
t610 = qJ(1,2) + pkin(7);
t597 = cos(t610);
t646 = 0.2e1 * qJ(3,2);
t613 = sin(t646);
t679 = pkin(7) + qJ(3,2);
t682 = -pkin(7) + qJ(3,2);
t711 = (t597 * t728 + sin(t610) * t729 + t629 * t730 + (sin(qJ(1,2) - t682) + sin(qJ(1,2) + t679)) * pkin(3)) / (t628 * t729 + pkin(3) * t613 + (sin(t679) + sin(t682)) * pkin(1));
t611 = qJ(1,1) + pkin(7);
t598 = cos(t611);
t647 = 0.2e1 * qJ(3,1);
t614 = sin(t647);
t680 = pkin(7) + qJ(3,1);
t683 = -pkin(7) + qJ(3,1);
t710 = (t598 * t728 + sin(t611) * t729 + t631 * t730 + (sin(qJ(1,1) - t683) + sin(qJ(1,1) + t680)) * pkin(3)) / (t630 * t729 + pkin(3) * t614 + (sin(t680) + sin(t683)) * pkin(1));
t623 = legFrame(3,2);
t585 = t623 + t609;
t586 = -t623 + t609;
t709 = -sin(t585) / 0.2e1 - sin(t586) / 0.2e1;
t624 = legFrame(2,2);
t587 = t624 + t610;
t588 = -t624 + t610;
t708 = -sin(t587) / 0.2e1 - sin(t588) / 0.2e1;
t625 = legFrame(1,2);
t589 = t625 + t611;
t590 = -t625 + t611;
t707 = -sin(t589) / 0.2e1 - sin(t590) / 0.2e1;
t706 = cos(t586) / 0.2e1 - cos(t585) / 0.2e1;
t705 = cos(t588) / 0.2e1 - cos(t587) / 0.2e1;
t704 = cos(t590) / 0.2e1 - cos(t589) / 0.2e1;
t570 = 0.1e1 / (t599 + t592);
t615 = 0.1e1 / t626;
t703 = t570 * t615;
t571 = 0.1e1 / (t599 + t593);
t616 = 0.1e1 / t628;
t702 = t571 * t616;
t572 = 0.1e1 / (t599 + t594);
t617 = 0.1e1 / t630;
t701 = t572 * t617;
t600 = sin(t623);
t697 = t626 * t600;
t603 = cos(t623);
t696 = t626 * t603;
t601 = sin(t624);
t693 = t628 * t601;
t604 = cos(t624);
t692 = t628 * t604;
t602 = sin(t625);
t689 = t630 * t602;
t605 = cos(t625);
t688 = t630 * t605;
t685 = t648 + t649;
t677 = pkin(3) * (-t633 * t622 + t695) * t632 ^ 2;
t676 = pkin(3) * (-t635 * t622 + t691) * t634 ^ 2;
t675 = pkin(3) * (-t637 * t622 + t687) * t636 ^ 2;
t674 = ((t592 * t627 + t633 * t642) * t622 - t583 * t627 + t621 * t700) * t615 * t632;
t673 = t600 * t715;
t672 = t603 * t715;
t671 = ((t593 * t629 + t635 * t642) * t622 - t583 * t629 + t621 * t699) * t616 * t634;
t670 = t601 * t714;
t669 = t604 * t714;
t668 = ((t594 * t631 + t637 * t642) * t622 - t583 * t631 + t621 * t698) * t617 * t636;
t667 = t602 * t713;
t666 = t605 * t713;
t665 = t650 * t712;
t664 = t650 * t711;
t663 = t650 * t710;
t638 = rSges(3,3) + pkin(5);
t662 = t638 + t722;
t661 = t715 * t725;
t660 = t714 * t724;
t659 = t713 * t723;
t560 = -t662 * t727 + Icges(3,6);
t561 = t662 * rSges(3,1) * m(3) - Icges(3,5);
t548 = t560 * t632 - t626 * t561;
t658 = t548 * t615 * t715;
t549 = t560 * t634 - t628 * t561;
t657 = t549 * t616 * t714;
t550 = t560 * t636 - t630 * t561;
t656 = t550 * t617 * t713;
t651 = pkin(1) ^ 2;
t652 = Icges(1,3) + Icges(2,3) + (0.2e1 * pkin(2) ^ 2 + 0.2e1 * t638 ^ 2 + 0.2e1 * t651 + t685) * t732 + t638 * t722 * t684 + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + t731 + Icges(3,1) / 0.2e1 + (t651 + (-0.2e1 * t722 + rSges(2,2)) * rSges(2,2) + (0.2e1 * t599 + rSges(2,1)) * rSges(2,1)) * m(2);
t643 = m(2) + m(3);
t595 = -rSges(3,1) * t727 + Icges(3,4);
t591 = t599 + pkin(2);
t582 = t685 * m(3) + Icges(3,3);
t538 = -t605 * t675 + (pkin(3) * t689 - t733 * t605) * t636 + t591 * t689;
t537 = -t604 * t676 + (pkin(3) * t693 - t734 * t604) * t634 + t591 * t693;
t536 = -t603 * t677 + (pkin(3) * t697 - t735 * t603) * t632 + t591 * t697;
t535 = t602 * t675 + (pkin(3) * t688 + t733 * t602) * t636 + t591 * t688;
t534 = t601 * t676 + (pkin(3) * t692 + t734 * t601) * t634 + t591 * t692;
t533 = t600 * t677 + (pkin(3) * t696 + t735 * t600) * t632 + t591 * t696;
t532 = cos(t647) * t726 + t595 * t614 + (t578 * t599 + (t599 + t578) * pkin(2)) * t684 + t652;
t531 = cos(t646) * t726 + t595 * t613 + (t577 * t599 + (t599 + t577) * pkin(2)) * t684 + t652;
t530 = cos(t645) * t726 + t595 * t612 + (t576 * t599 + (t599 + t576) * pkin(2)) * t684 + t652;
t529 = -t572 * t643 * t668 + t663 * t723;
t528 = -t571 * t643 * t671 + t664 * t724;
t527 = -t570 * t643 * t674 + t665 * t725;
t526 = (t538 * t643 - t605 * t659) * t701;
t525 = (t537 * t643 - t604 * t660) * t702;
t524 = (t536 * t643 - t603 * t661) * t703;
t523 = (t535 * t643 + t602 * t659) * t701;
t522 = (t534 * t643 + t601 * t660) * t702;
t521 = (t533 * t643 + t600 * t661) * t703;
t520 = t582 * t663 + (-t550 * t598 - t668 * t723) * t572;
t519 = t582 * t664 + (-t549 * t597 - t671 * t724) * t571;
t518 = t582 * t665 + (-t548 * t596 - t674 * t725) * t570;
t517 = -t598 * t572 * t532 + t550 * t663;
t516 = -t597 * t571 * t531 + t549 * t664;
t515 = -t596 * t570 * t530 + t548 * t665;
t514 = (t532 * t704 + t602 * t656) * t572;
t513 = (t531 * t705 + t601 * t657) * t571;
t512 = (t530 * t706 + t600 * t658) * t570;
t511 = (t532 * t707 - t605 * t656) * t572;
t510 = (t531 * t708 - t604 * t657) * t571;
t509 = (t530 * t709 - t603 * t658) * t570;
t508 = (t550 * t704 + (t535 * t723 + t582 * t667) * t617) * t572;
t507 = (t549 * t705 + (t534 * t724 + t582 * t670) * t616) * t571;
t506 = (t548 * t706 + (t533 * t725 + t582 * t673) * t615) * t570;
t505 = (t550 * t707 + (t538 * t723 - t582 * t666) * t617) * t572;
t504 = (t549 * t708 + (t537 * t724 - t582 * t669) * t616) * t571;
t503 = (t548 * t709 + (t536 * t725 - t582 * t672) * t615) * t570;
t1 = [m(4) + (t511 * t707 + (-t505 * t666 + t526 * t538) * t617) * t572 + (t510 * t708 + (-t504 * t669 + t525 * t537) * t616) * t571 + (t509 * t709 + (-t503 * t672 + t524 * t536) * t615) * t570, (t511 * t704 + (t505 * t667 + t526 * t535) * t617) * t572 + (t510 * t705 + (t504 * t670 + t525 * t534) * t616) * t571 + (t509 * t706 + (t503 * t673 + t524 * t533) * t615) * t570, (-t511 * t598 - t526 * t668) * t572 + (-t510 * t597 - t525 * t671) * t571 + (-t509 * t596 - t524 * t674) * t570 + (t503 * t712 + t504 * t711 + t505 * t710) * t650; (t514 * t707 + (-t508 * t666 + t523 * t538) * t617) * t572 + (t513 * t708 + (-t507 * t669 + t522 * t537) * t616) * t571 + (t512 * t709 + (-t506 * t672 + t521 * t536) * t615) * t570, m(4) + (t514 * t704 + (t508 * t667 + t523 * t535) * t617) * t572 + (t513 * t705 + (t507 * t670 + t522 * t534) * t616) * t571 + (t512 * t706 + (t506 * t673 + t521 * t533) * t615) * t570, (-t514 * t598 - t523 * t668) * t572 + (-t513 * t597 - t522 * t671) * t571 + (-t512 * t596 - t521 * t674) * t570 + (t506 * t712 + t507 * t711 + t508 * t710) * t650; (t517 * t707 + (-t520 * t666 + t529 * t538) * t617) * t572 + (t516 * t708 + (-t519 * t669 + t528 * t537) * t616) * t571 + (t515 * t709 + (-t518 * t672 + t527 * t536) * t615) * t570, (t517 * t704 + (t520 * t667 + t529 * t535) * t617) * t572 + (t516 * t705 + (t519 * t670 + t528 * t534) * t616) * t571 + (t515 * t706 + (t518 * t673 + t527 * t533) * t615) * t570, m(4) + (-t517 * t598 - t529 * t668) * t572 + (-t516 * t597 - t528 * t671) * t571 + (-t515 * t596 - t527 * t674) * t570 + (t518 * t712 + t519 * t711 + t520 * t710) * t650;];
MX  = t1;
