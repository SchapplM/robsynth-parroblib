% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:10:32
% EndTime: 2020-03-02 20:10:39
% DurationCPUTime: 7.45s
% Computational Cost: add. (23509->267), mult. (15450->521), div. (3248->6), fcn. (15916->26), ass. (0->260)
t726 = 0.1e1 / pkin(2);
t701 = pkin(7) + qJ(2,1);
t685 = qJ(3,1) + t701;
t669 = sin(t685);
t672 = cos(t685);
t679 = sin(t701);
t682 = cos(t701);
t821 = 0.1e1 / (t669 * t682 - t679 * t672);
t837 = t821 * t726;
t700 = pkin(7) + qJ(2,2);
t684 = qJ(3,2) + t700;
t668 = sin(t684);
t671 = cos(t684);
t678 = sin(t700);
t681 = cos(t700);
t822 = 0.1e1 / (t668 * t681 - t678 * t671);
t836 = t822 * t726;
t699 = pkin(7) + qJ(2,3);
t683 = qJ(3,3) + t699;
t667 = sin(t683);
t670 = cos(t683);
t677 = sin(t699);
t680 = cos(t699);
t823 = 0.1e1 / (t667 * t680 - t677 * t670);
t835 = t823 * t726;
t698 = pkin(7) + qJ(2,4);
t676 = qJ(3,4) + t698;
t665 = sin(t676);
t666 = cos(t676);
t674 = sin(t698);
t675 = cos(t698);
t824 = 0.1e1 / (t665 * t675 - t674 * t666);
t834 = t824 * t726;
t706 = legFrame(1,3);
t689 = sin(t706);
t693 = cos(t706);
t626 = t693 * t669 + t689 * t672;
t627 = -t669 * t689 + t693 * t672;
t710 = xP(4);
t696 = sin(t710);
t697 = cos(t710);
t716 = koppelP(1,2);
t720 = koppelP(1,1);
t651 = -t696 * t716 + t697 * t720;
t707 = xDP(4);
t708 = xDP(2);
t631 = t651 * t707 + t708;
t647 = t696 * t720 + t697 * t716;
t709 = xDP(1);
t635 = -t647 * t707 + t709;
t571 = (t626 * t631 + t627 * t635) * t837;
t724 = 0.1e1 / pkin(3);
t786 = t724 * t726;
t764 = t821 * t786;
t599 = -pkin(2) * (t679 * t689 - t693 * t682) + t627 * pkin(3);
t808 = t599 * t635;
t596 = pkin(2) * (t693 * t679 + t689 * t682) + t626 * pkin(3);
t811 = t596 * t631;
t536 = -(t808 / 0.2e1 + t811 / 0.2e1) * t764 + t571;
t727 = (t808 + t811) * t764;
t539 = -t727 + t571;
t723 = pkin(3) ^ 2;
t725 = pkin(2) ^ 2;
t743 = t669 * t679 + t672 * t682;
t731 = t743 * pkin(2);
t780 = 0.2e1 * pkin(2) * pkin(3);
t816 = t539 * t727;
t519 = ((-t743 * t536 * t780 - t539 * t723 - t725 * t571) * t724 * t571 + (pkin(3) + t731) * t816) * t837;
t523 = (-pkin(3) * t816 + (pkin(3) * t539 + t571 * t731) * t571) * t837;
t694 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t642 = -rSges(3,1) * t682 - t679 * rSges(3,2);
t643 = t679 * rSges(3,1) - rSges(3,2) * t682;
t735 = (-t642 * t672 + t643 * t669) * pkin(2);
t591 = Icges(3,3) + (t694 + t735) * m(3);
t652 = t694 * m(3) + Icges(3,3);
t515 = t652 * t519 + t591 * t523;
t604 = t642 * t669 + t643 * t672;
t781 = m(3) * t571 ^ 2 * t604;
t789 = t821 * t724;
t833 = -t515 * t764 + t781 * t789;
t705 = legFrame(2,3);
t688 = sin(t705);
t692 = cos(t705);
t624 = t692 * t668 + t688 * t671;
t625 = -t668 * t688 + t692 * t671;
t715 = koppelP(2,2);
t719 = koppelP(2,1);
t650 = -t696 * t715 + t697 * t719;
t630 = t650 * t707 + t708;
t646 = t696 * t719 + t697 * t715;
t634 = -t646 * t707 + t709;
t570 = (t624 * t630 + t625 * t634) * t836;
t766 = t822 * t786;
t598 = -pkin(2) * (t678 * t688 - t692 * t681) + t625 * pkin(3);
t809 = t598 * t634;
t595 = pkin(2) * (t692 * t678 + t688 * t681) + t624 * pkin(3);
t812 = t595 * t630;
t535 = -(t809 / 0.2e1 + t812 / 0.2e1) * t766 + t570;
t728 = (t809 + t812) * t766;
t538 = -t728 + t570;
t744 = t668 * t678 + t671 * t681;
t732 = t744 * pkin(2);
t817 = t538 * t728;
t518 = ((-t744 * t535 * t780 - t538 * t723 - t725 * t570) * t724 * t570 + (pkin(3) + t732) * t817) * t836;
t522 = (-pkin(3) * t817 + (pkin(3) * t538 + t570 * t732) * t570) * t836;
t640 = -rSges(3,1) * t681 - t678 * rSges(3,2);
t641 = t678 * rSges(3,1) - rSges(3,2) * t681;
t736 = (-t640 * t671 + t641 * t668) * pkin(2);
t590 = Icges(3,3) + (t694 + t736) * m(3);
t514 = t652 * t518 + t590 * t522;
t603 = t640 * t668 + t641 * t671;
t782 = m(3) * t570 ^ 2 * t603;
t791 = t822 * t724;
t832 = -t514 * t766 + t782 * t791;
t704 = legFrame(3,3);
t687 = sin(t704);
t691 = cos(t704);
t622 = t691 * t667 + t687 * t670;
t623 = -t667 * t687 + t691 * t670;
t714 = koppelP(3,2);
t718 = koppelP(3,1);
t649 = -t696 * t714 + t697 * t718;
t629 = t649 * t707 + t708;
t645 = t696 * t718 + t697 * t714;
t633 = -t645 * t707 + t709;
t569 = (t622 * t629 + t623 * t633) * t835;
t768 = t823 * t786;
t597 = -pkin(2) * (t677 * t687 - t691 * t680) + t623 * pkin(3);
t810 = t597 * t633;
t594 = pkin(2) * (t691 * t677 + t687 * t680) + t622 * pkin(3);
t813 = t594 * t629;
t534 = -(t810 / 0.2e1 + t813 / 0.2e1) * t768 + t569;
t729 = (t810 + t813) * t768;
t537 = -t729 + t569;
t745 = t667 * t677 + t670 * t680;
t733 = t745 * pkin(2);
t818 = t537 * t729;
t517 = ((-t745 * t534 * t780 - t537 * t723 - t725 * t569) * t724 * t569 + (pkin(3) + t733) * t818) * t835;
t521 = (-pkin(3) * t818 + (pkin(3) * t537 + t569 * t733) * t569) * t835;
t638 = -rSges(3,1) * t680 - t677 * rSges(3,2);
t639 = t677 * rSges(3,1) - rSges(3,2) * t680;
t737 = (-t638 * t670 + t639 * t667) * pkin(2);
t589 = Icges(3,3) + (t694 + t737) * m(3);
t513 = t652 * t517 + t589 * t521;
t602 = t638 * t667 + t639 * t670;
t783 = m(3) * t569 ^ 2 * t602;
t793 = t823 * t724;
t831 = -t513 * t768 + t783 * t793;
t703 = legFrame(4,3);
t686 = sin(t703);
t690 = cos(t703);
t620 = t690 * t665 + t686 * t666;
t621 = -t665 * t686 + t690 * t666;
t713 = koppelP(4,2);
t717 = koppelP(4,1);
t648 = -t696 * t713 + t697 * t717;
t628 = t648 * t707 + t708;
t644 = t696 * t717 + t697 * t713;
t632 = -t644 * t707 + t709;
t565 = (t620 * t628 + t621 * t632) * t834;
t770 = t824 * t786;
t593 = -pkin(2) * (t674 * t686 - t690 * t675) + t621 * pkin(3);
t814 = t593 * t632;
t592 = pkin(2) * (t690 * t674 + t686 * t675) + t620 * pkin(3);
t815 = t592 * t628;
t532 = -(t814 / 0.2e1 + t815 / 0.2e1) * t770 + t565;
t730 = (t814 + t815) * t770;
t533 = -t730 + t565;
t746 = t665 * t674 + t666 * t675;
t734 = t746 * pkin(2);
t819 = t533 * t730;
t516 = ((-t746 * t532 * t780 - t533 * t723 - t725 * t565) * t724 * t565 + (pkin(3) + t734) * t819) * t834;
t520 = (-pkin(3) * t819 + (pkin(3) * t533 + t565 * t734) * t565) * t834;
t636 = -rSges(3,1) * t675 - t674 * rSges(3,2);
t637 = t674 * rSges(3,1) - rSges(3,2) * t675;
t738 = (-t636 * t666 + t637 * t665) * pkin(2);
t588 = Icges(3,3) + (t694 + t738) * m(3);
t509 = t652 * t516 + t588 * t520;
t600 = t636 * t665 + t637 * t666;
t784 = m(3) * t565 ^ 2 * t600;
t804 = t824 * t724;
t830 = -t509 * t770 + t784 * t804;
t673 = t725 + t694;
t751 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3) + Icges(3,3);
t587 = (t673 + 0.2e1 * t735) * m(3) + t751;
t512 = t591 * t519 + t587 * t523;
t825 = -0.2e1 * m(3);
t752 = t536 * t727 * t604 * t825;
t829 = t512 * t837 + t821 * t752;
t586 = (t673 + 0.2e1 * t736) * m(3) + t751;
t511 = t590 * t518 + t586 * t522;
t753 = t535 * t728 * t603 * t825;
t828 = t511 * t836 + t822 * t753;
t585 = (t673 + 0.2e1 * t737) * m(3) + t751;
t510 = t589 * t517 + t585 * t521;
t754 = t534 * t729 * t602 * t825;
t827 = t510 * t835 + t823 * t754;
t584 = (t673 + 0.2e1 * t738) * m(3) + t751;
t508 = t588 * t516 + t584 * t520;
t755 = t532 * t730 * t600 * t825;
t826 = t508 * t834 + t824 * t755;
t702 = t707 ^ 2;
t820 = m(4) * t702;
t807 = t824 * t620;
t806 = t824 * t621;
t802 = t823 * t622;
t801 = t823 * t623;
t799 = t822 * t624;
t798 = t822 * t625;
t796 = t821 * t626;
t795 = t821 * t627;
t787 = t652 * t724;
t785 = t726 * t702;
t779 = t592 * t804;
t778 = t593 * t804;
t777 = t594 * t793;
t776 = t595 * t791;
t775 = t596 * t789;
t774 = t597 * t793;
t773 = t598 * t791;
t772 = t599 * t789;
t771 = t824 * t787;
t769 = t823 * t787;
t767 = t822 * t787;
t765 = t821 * t787;
t712 = rSges(4,1);
t711 = rSges(4,2);
t575 = (t626 * t651 - t627 * t647) * t837;
t574 = (t624 * t650 - t625 * t646) * t836;
t573 = (t622 * t649 - t623 * t645) * t835;
t572 = (t620 * t648 - t621 * t644) * t834;
t563 = (t596 * t651 - t599 * t647) * t764;
t562 = (t595 * t650 - t598 * t646) * t766;
t561 = (t594 * t649 - t597 * t645) * t768;
t560 = (t592 * t648 - t593 * t644) * t770;
t559 = (t591 * t795 - t599 * t765) * t726;
t558 = (t591 * t796 - t596 * t765) * t726;
t557 = (t590 * t798 - t598 * t767) * t726;
t556 = (t590 * t799 - t595 * t767) * t726;
t555 = (t589 * t801 - t597 * t769) * t726;
t554 = (t589 * t802 - t594 * t769) * t726;
t550 = (t588 * t806 - t593 * t771) * t726;
t549 = (t588 * t807 - t592 * t771) * t726;
t547 = (t587 * t795 - t591 * t772) * t726;
t546 = (t587 * t796 - t591 * t775) * t726;
t545 = (t586 * t798 - t590 * t773) * t726;
t544 = (t586 * t799 - t590 * t776) * t726;
t543 = (t585 * t801 - t589 * t774) * t726;
t542 = (t585 * t802 - t589 * t777) * t726;
t541 = (t584 * t806 - t588 * t778) * t726;
t540 = (t584 * t807 - t588 * t779) * t726;
t531 = -t563 * t652 + t575 * t591;
t530 = -t562 * t652 + t574 * t590;
t529 = -t561 * t652 + t573 * t589;
t528 = -t560 * t652 + t572 * t588;
t527 = -t563 * t591 + t575 * t587;
t526 = -t562 * t590 + t574 * t586;
t525 = -t561 * t589 + t573 * t585;
t524 = -t560 * t588 + t572 * t584;
t1 = [(t696 * t711 - t697 * t712) * t820 + t829 * t627 + t828 * t625 + t827 * t623 + t826 * t621 + t833 * t599 + t832 * t598 + t831 * t597 + t830 * t593 + (-(t547 * t795 - t559 * t772) * t651 - (t547 * t796 - t559 * t775) * t647 - (t545 * t798 - t557 * t773) * t650 - (t545 * t799 - t557 * t776) * t646 - (t543 * t801 - t555 * t774) * t649 - (t543 * t802 - t555 * t777) * t645 - (t541 * t806 - t550 * t778) * t648 - (t541 * t807 - t550 * t779) * t644) * t785; -(t696 * t712 + t697 * t711) * t820 + t829 * t626 + t828 * t624 + t827 * t622 + t826 * t620 + t833 * t596 + t832 * t595 + t831 * t594 + t830 * t592 + (-(t546 * t795 - t558 * t772) * t651 - (t546 * t796 - t558 * t775) * t647 - (t544 * t798 - t556 * t773) * t650 - (t544 * t799 - t556 * t776) * t646 - (t542 * t801 - t554 * t774) * t649 - (t542 * t802 - t554 * t777) * t645 - (t540 * t806 - t549 * t778) * t648 - (t540 * t807 - t549 * t779) * t644) * t785; 0; t572 * t508 - t560 * t509 + t573 * t510 + t574 * t511 + t575 * t512 - t561 * t513 - t562 * t514 - t563 * t515 + (-(t527 * t795 - t531 * t772) * t651 - (t527 * t796 - t531 * t775) * t647 - (t526 * t798 - t530 * t773) * t650 - (t526 * t799 - t530 * t776) * t646 - (t525 * t801 - t529 * t774) * t649 - (t525 * t802 - t529 * t777) * t645 - (t524 * t806 - t528 * t778) * t648 - (t524 * t807 - t528 * t779) * t644) * t785 + (t560 * t784 + t561 * t783 + t562 * t782 + t563 * t781 + t572 * t755 + t573 * t754 + t574 * t753 + t575 * t752) * pkin(2);];
taucX  = t1;
