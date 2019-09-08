% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPP1G1P1A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:01
% EndTime: 2019-05-03 14:52:06
% DurationCPUTime: 4.60s
% Computational Cost: add. (27568->346), mult. (37633->547), div. (2085->3), fcn. (14860->14), ass. (0->265)
t658 = (qJ(3,3) ^ 2);
t675 = (pkin(1) ^ 2);
t753 = 1 + t675;
t755 = 2 * qJ(3,3);
t609 = pkin(1) * t755 + t658 + t753;
t645 = sin(qJ(1,3));
t642 = pkin(1) + qJ(3,3);
t648 = cos(qJ(1,3));
t717 = t642 * t648;
t704 = qJ(2,3) * t717;
t579 = t609 * t645 - t704;
t718 = t642 * t645;
t612 = qJ(2,3) * t718;
t582 = t609 * t648 + t612;
t639 = legFrame(3,3);
t615 = sin(t639);
t618 = cos(t639);
t552 = t579 * t618 + t582 * t615;
t555 = -t579 * t615 + t582 * t618;
t655 = xP(3);
t623 = sin(t655);
t624 = cos(t655);
t664 = koppelP(3,2);
t667 = koppelP(3,1);
t594 = -t623 * t664 + t624 * t667;
t651 = xDP(3);
t652 = xDP(2);
t573 = t594 * t651 + t652;
t591 = t623 * t667 + t624 * t664;
t653 = xDP(1);
t576 = -t591 * t651 + t653;
t659 = qJ(2,3) ^ 2;
t768 = t609 + t659;
t606 = 0.1e1 / t768;
t516 = (t552 * t573 + t555 * t576) * t606;
t777 = t516 * t642;
t660 = (qJ(3,2) ^ 2);
t757 = 2 * qJ(3,2);
t610 = pkin(1) * t757 + t660 + t753;
t646 = sin(qJ(1,2));
t643 = pkin(1) + qJ(3,2);
t649 = cos(qJ(1,2));
t715 = t643 * t649;
t705 = qJ(2,2) * t715;
t580 = t610 * t646 - t705;
t716 = t643 * t646;
t613 = qJ(2,2) * t716;
t583 = t610 * t649 + t613;
t640 = legFrame(2,3);
t616 = sin(t640);
t619 = cos(t640);
t553 = t580 * t619 + t583 * t616;
t556 = -t580 * t616 + t583 * t619;
t665 = koppelP(2,2);
t668 = koppelP(2,1);
t595 = -t623 * t665 + t624 * t668;
t574 = t595 * t651 + t652;
t592 = t623 * t668 + t624 * t665;
t577 = -t592 * t651 + t653;
t661 = qJ(2,2) ^ 2;
t767 = t610 + t661;
t607 = 0.1e1 / t767;
t517 = (t553 * t574 + t556 * t577) * t607;
t776 = t517 * t643;
t662 = (qJ(3,1) ^ 2);
t759 = 2 * qJ(3,1);
t611 = pkin(1) * t759 + t662 + t753;
t647 = sin(qJ(1,1));
t644 = pkin(1) + qJ(3,1);
t650 = cos(qJ(1,1));
t713 = t644 * t650;
t706 = qJ(2,1) * t713;
t581 = t611 * t647 - t706;
t714 = t644 * t647;
t614 = qJ(2,1) * t714;
t584 = t611 * t650 + t614;
t641 = legFrame(1,3);
t617 = sin(t641);
t620 = cos(t641);
t554 = t581 * t620 + t584 * t617;
t557 = -t581 * t617 + t584 * t620;
t666 = koppelP(1,2);
t669 = koppelP(1,1);
t596 = -t623 * t666 + t624 * t669;
t575 = t596 * t651 + t652;
t593 = t623 * t669 + t624 * t666;
t578 = -t593 * t651 + t653;
t663 = qJ(2,1) ^ 2;
t766 = t611 + t663;
t608 = 0.1e1 / t766;
t518 = (t554 * t575 + t557 * t578) * t608;
t775 = t518 * t644;
t628 = 0.1e1 + t659;
t629 = 0.1e1 + t661;
t630 = 0.1e1 + t663;
t599 = qJ(2,1) * t650 - t714;
t602 = qJ(2,1) * t647 + t713;
t569 = t599 * t620 - t602 * t617;
t572 = t599 * t617 + t602 * t620;
t771 = (t569 * t578 + t572 * t575) * t608;
t598 = qJ(2,2) * t649 - t716;
t601 = qJ(2,2) * t646 + t715;
t568 = t598 * t619 - t601 * t616;
t571 = t598 * t616 + t601 * t619;
t770 = (t568 * t577 + t571 * t574) * t607;
t597 = qJ(2,3) * t648 - t718;
t600 = qJ(2,3) * t645 + t717;
t567 = t597 * t618 - t600 * t615;
t570 = t597 * t615 + t600 * t618;
t769 = (t567 * t576 + t570 * t573) * t606;
t765 = -t658 - t659;
t764 = -t660 - t661;
t763 = -t662 - t663;
t762 = 2 * m(3);
t761 = 2 * pkin(1);
t632 = t651 ^ 2;
t760 = 0.2e1 * qJ(2,1);
t758 = 0.2e1 * qJ(2,2);
t756 = 0.2e1 * qJ(2,3);
t754 = -3 * t675;
t752 = m(2) * (rSges(2,3) + qJ(2,3));
t751 = m(2) * (rSges(2,3) + qJ(2,2));
t750 = m(2) * (rSges(2,3) + qJ(2,1));
t749 = m(3) * t606;
t748 = m(3) * t607;
t747 = m(3) * t608;
t636 = rSges(3,2) + qJ(2,3);
t746 = m(3) * t636;
t637 = rSges(3,2) + qJ(2,2);
t745 = m(3) * t637;
t638 = rSges(3,2) + qJ(2,1);
t744 = m(3) * t638;
t743 = m(4) * t632;
t742 = rSges(3,3) + qJ(3,1);
t741 = rSges(3,3) + qJ(3,2);
t740 = rSges(3,3) + qJ(3,3);
t587 = t630 * t647 + t706;
t590 = -t630 * t650 + t614;
t563 = t587 * t620 - t590 * t617;
t548 = t563 * t608 * t578;
t566 = t587 * t617 + t590 * t620;
t551 = t566 * t608 * t575;
t524 = t548 + t551;
t738 = qJ(2,1) * t524;
t586 = t629 * t646 + t705;
t589 = -t629 * t649 + t613;
t562 = t586 * t619 - t589 * t616;
t547 = t562 * t607 * t577;
t565 = t586 * t616 + t589 * t619;
t550 = t565 * t607 * t574;
t523 = t547 + t550;
t736 = qJ(2,2) * t523;
t585 = t628 * t645 + t704;
t588 = -t628 * t648 + t612;
t561 = t585 * t618 - t588 * t615;
t546 = t561 * t606 * t576;
t564 = t585 * t615 + t588 * t618;
t549 = t564 * t606 * t573;
t522 = t546 + t549;
t734 = qJ(2,3) * t522;
t733 = t769 * t606;
t732 = t770 * t607;
t731 = t771 * t608;
t674 = pkin(1) * t675;
t480 = (0.2e1 * t609 * t522 - 0.2e1 * qJ(2,3) * t777 + (-t674 + (t754 + t765) * qJ(3,3) - qJ(3,3) + (-(3 * t658) - t628) * pkin(1)) * t769) * t733;
t486 = (0.2e1 * t642 * t734 - t516 + (-t659 - t628) * t516 - t768 * t769 * qJ(2,3)) * t733;
t489 = 0.2e1 * (t734 + t777) * t733;
t686 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,1) + Icges(1,3);
t693 = rSges(2,3) ^ 2 + t675 + ((-2 * pkin(1) + rSges(2,2)) * rSges(2,2));
t700 = rSges(3,2) ^ 2 + (rSges(3,3) ^ 2) + t675;
t558 = (rSges(3,2) * t756 + (rSges(3,3) * t755) + (t740 * t761) + t700 - t765) * m(3) + (rSges(2,3) * t756 + t659 + t693) * m(2) + t686;
t622 = (pkin(1) - rSges(2,2)) * m(2);
t625 = -pkin(1) - t740;
t603 = (m(3) * t625) - t622;
t474 = -t480 * t746 - t486 * t603 - t489 * t558;
t730 = t606 * t474;
t477 = (-t489 * t636 - t480) * m(3);
t729 = t606 * t477;
t654 = m(2) + m(3);
t483 = -t486 * t654 - t489 * t603;
t728 = t606 * t483;
t727 = t606 * t632;
t481 = (0.2e1 * t610 * t523 - 0.2e1 * qJ(2,2) * t776 + (-t674 + (t754 + t764) * qJ(3,2) - qJ(3,2) + (-(3 * t660) - t629) * pkin(1)) * t770) * t732;
t487 = (0.2e1 * t643 * t736 - t517 + (-t661 - t629) * t517 - t767 * t770 * qJ(2,2)) * t732;
t490 = 0.2e1 * (t736 + t776) * t732;
t559 = (rSges(3,2) * t758 + (rSges(3,3) * t757) + (t741 * t761) + t700 - t764) * m(3) + (rSges(2,3) * t758 + t661 + t693) * m(2) + t686;
t626 = -pkin(1) - t741;
t604 = (m(3) * t626) - t622;
t475 = -t481 * t745 - t487 * t604 - t490 * t559;
t726 = t607 * t475;
t478 = (-t490 * t637 - t481) * m(3);
t725 = t607 * t478;
t484 = -t487 * t654 - t490 * t604;
t724 = t607 * t484;
t723 = t607 * t632;
t482 = (0.2e1 * t611 * t524 - 0.2e1 * qJ(2,1) * t775 + (-t674 + (t754 + t763) * qJ(3,1) - qJ(3,1) + (-(3 * t662) - t630) * pkin(1)) * t771) * t731;
t488 = (0.2e1 * t644 * t738 - t518 + (-t663 - t630) * t518 - t766 * t771 * qJ(2,1)) * t731;
t491 = 0.2e1 * (t738 + t775) * t731;
t560 = (rSges(3,2) * t760 + (rSges(3,3) * t759) + (t742 * t761) + t700 - t763) * m(3) + (rSges(2,3) * t760 + t663 + t693) * m(2) + t686;
t627 = -pkin(1) - t742;
t605 = (m(3) * t627) - t622;
t476 = -t482 * t744 - t488 * t605 - t491 * t560;
t722 = t608 * t476;
t479 = (-t491 * t638 - t482) * m(3);
t721 = t608 * t479;
t485 = -t488 * t654 - t491 * t605;
t720 = t608 * t485;
t719 = t608 * t632;
t712 = 0.2e1 * ((-t516 * t625 + t522 * t636) * m(3) + t522 * t752) * t769;
t711 = 0.2e1 * ((-t517 * t626 + t523 * t637) * m(3) + t523 * t751) * t770;
t710 = 0.2e1 * ((-t518 * t627 + t524 * t638) * m(3) + t524 * t750) * t771;
t709 = m(3) * (-t625 * t769 - 0.2e1 * t546 - 0.2e1 * t549) * t769;
t708 = m(3) * (-t626 * t770 - 0.2e1 * t547 - 0.2e1 * t550) * t770;
t707 = m(3) * (-t627 * t771 - 0.2e1 * t548 - 0.2e1 * t551) * t771;
t498 = (t746 + t752) * t769 + t516 * t762;
t703 = t498 * t733;
t499 = (t745 + t751) * t770 + t517 * t762;
t702 = t499 * t732;
t500 = (t744 + t750) * t771 + t518 * t762;
t701 = t500 * t731;
t699 = t606 * t712;
t698 = t607 * t711;
t697 = t608 * t710;
t696 = t606 * t709;
t695 = t607 * t708;
t694 = t608 * t707;
t657 = rSges(4,1);
t656 = rSges(4,2);
t545 = (t566 * t654 + t572 * t605) * t608;
t544 = (t565 * t654 + t571 * t604) * t607;
t543 = (t564 * t654 + t570 * t603) * t606;
t542 = (t563 * t654 + t569 * t605) * t608;
t541 = (t562 * t654 + t568 * t604) * t607;
t540 = (t561 * t654 + t567 * t603) * t606;
t539 = (t572 * t638 + t554) * t747;
t538 = (t571 * t637 + t553) * t748;
t537 = (t570 * t636 + t552) * t749;
t536 = (t569 * t638 + t557) * t747;
t535 = (t568 * t637 + t556) * t748;
t534 = (t567 * t636 + t555) * t749;
t533 = (-t569 * t593 + t572 * t596) * t608;
t532 = (-t568 * t592 + t571 * t595) * t607;
t531 = (-t567 * t591 + t570 * t594) * t606;
t527 = (-t563 * t593 + t566 * t596) * t608;
t526 = (-t562 * t592 + t565 * t595) * t607;
t525 = (-t561 * t591 + t564 * t594) * t606;
t521 = (t554 * t596 - t557 * t593) * t608;
t520 = (t553 * t595 - t556 * t592) * t607;
t519 = (t552 * t594 - t555 * t591) * t606;
t515 = (t554 * t744 + t560 * t572 + t566 * t605) * t608;
t514 = (t553 * t745 + t559 * t571 + t565 * t604) * t607;
t513 = (t552 * t746 + t558 * t570 + t564 * t603) * t606;
t512 = (t557 * t744 + t560 * t569 + t563 * t605) * t608;
t511 = (t556 * t745 + t559 * t568 + t562 * t604) * t607;
t510 = (t555 * t746 + t558 * t567 + t561 * t603) * t606;
t509 = t527 * t654 + t533 * t605;
t508 = t526 * t654 + t532 * t604;
t507 = t525 * t654 + t531 * t603;
t506 = (t533 * t638 + t521) * m(3);
t505 = (t532 * t637 + t520) * m(3);
t504 = (t531 * t636 + t519) * m(3);
t497 = t521 * t744 + t527 * t605 + t533 * t560;
t496 = t520 * t745 + t526 * t604 + t532 * t559;
t495 = t519 * t746 + t525 * t603 + t531 * t558;
t1 = [(-(t512 * t569 + t536 * t557 + t542 * t563) * t596 - (t512 * t572 + t536 * t554 + t542 * t566) * t593) * t719 + t569 * t722 + t563 * t720 + t557 * t721 + t569 * t697 - t563 * t701 - t557 * t694 + (-(t511 * t568 + t535 * t556 + t541 * t562) * t595 - (t511 * t571 + t535 * t553 + t541 * t565) * t592) * t723 + t568 * t726 + t562 * t724 + t556 * t725 + t568 * t698 - t562 * t702 - t556 * t695 + (-(t510 * t567 + t534 * t555 + t540 * t561) * t594 - (t510 * t570 + t534 * t552 + t540 * t564) * t591) * t727 + t567 * t730 + t561 * t728 + t555 * t729 + t567 * t699 - t561 * t703 - t555 * t696 - (-t623 * t656 + t624 * t657) * t743; (-(t515 * t569 + t539 * t557 + t545 * t563) * t596 - (t515 * t572 + t539 * t554 + t545 * t566) * t593) * t719 + t572 * t722 + t566 * t720 + t554 * t721 + t572 * t697 - t566 * t701 - t554 * t694 + (-(t514 * t568 + t538 * t556 + t544 * t562) * t595 - (t514 * t571 + t538 * t553 + t544 * t565) * t592) * t723 + t571 * t726 + t565 * t724 + t553 * t725 + t571 * t698 - t565 * t702 - t553 * t695 + (-(t513 * t567 + t537 * t555 + t543 * t561) * t594 - (t513 * t570 + t537 * t552 + t543 * t564) * t591) * t727 + t570 * t730 + t564 * t728 + t552 * t729 + t570 * t699 - t564 * t703 - t552 * t696 - (t623 * t657 + t624 * t656) * t743; (-(t497 * t569 + t506 * t557 + t509 * t563) * t596 - (t497 * t572 + t506 * t554 + t509 * t566) * t593) * t719 + (-(t496 * t568 + t505 * t556 + t508 * t562) * t595 - (t496 * t571 + t505 * t553 + t508 * t565) * t592) * t723 + (-(t495 * t567 + t504 * t555 + t507 * t561) * t594 - (t495 * t570 + t504 * t552 + t507 * t564) * t591) * t727 + (t476 + t710) * t533 + (t475 + t711) * t532 + (t474 + t712) * t531 + (-t500 * t771 + t485) * t527 + (-t499 * t770 + t484) * t526 + (-t498 * t769 + t483) * t525 + (t479 - t707) * t521 + (t478 - t708) * t520 + (t477 - t709) * t519;];
taucX  = t1;
