% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G4A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 18:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G4A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:14:21
% EndTime: 2020-08-06 18:14:24
% DurationCPUTime: 2.61s
% Computational Cost: add. (9354->343), mult. (18738->666), div. (432->10), fcn. (19692->34), ass. (0->311)
t827 = (m(3) * rSges(3,2));
t828 = -2 * rSges(3,1) * t827 + 2 * Icges(3,4);
t703 = (pkin(2) * m(3));
t678 = cos(pkin(4));
t688 = sin(qJ(3,3));
t694 = cos(qJ(3,3));
t715 = rSges(3,1) * t694 - rSges(3,2) * t688;
t676 = sin(pkin(4));
t689 = sin(qJ(2,3));
t757 = t676 * t689;
t826 = m(3) * (t715 * t678 - (t688 * rSges(3,1) + t694 * rSges(3,2)) * t757);
t690 = sin(qJ(3,2));
t696 = cos(qJ(3,2));
t714 = rSges(3,1) * t696 - rSges(3,2) * t690;
t691 = sin(qJ(2,2));
t756 = t676 * t691;
t825 = m(3) * (t714 * t678 - (t690 * rSges(3,1) + t696 * rSges(3,2)) * t756);
t692 = sin(qJ(3,1));
t698 = cos(qJ(3,1));
t713 = rSges(3,1) * t698 - rSges(3,2) * t692;
t693 = sin(qJ(2,1));
t755 = t676 * t693;
t824 = m(3) * (t713 * t678 - (t692 * rSges(3,1) + t698 * rSges(3,2)) * t755);
t672 = t694 ^ 2;
t823 = pkin(3) * t672;
t673 = t696 ^ 2;
t822 = pkin(3) * t673;
t674 = t698 ^ 2;
t821 = pkin(3) * t674;
t820 = pkin(3) * t676;
t819 = t688 * pkin(2);
t818 = t690 * pkin(2);
t817 = t692 * pkin(2);
t816 = t694 * pkin(3);
t815 = t696 * pkin(3);
t814 = t698 * pkin(3);
t700 = pkin(6) + rSges(3,3);
t813 = t700 * m(3);
t695 = cos(qJ(2,3));
t702 = pkin(7) + pkin(6);
t747 = t689 * t702;
t632 = pkin(2) * t695 + t747;
t675 = sin(pkin(8));
t677 = cos(pkin(8));
t645 = t702 * t695;
t629 = pkin(2) * t689 - t645;
t709 = -t629 * t678 + t688 * t820;
t581 = -t632 * t677 - t709 * t675;
t584 = t632 * t675 - t709 * t677;
t748 = t688 * t678;
t602 = pkin(3) * t748 + t629 * t676;
t679 = legFrame(3,3);
t649 = sin(t679);
t655 = cos(t679);
t685 = legFrame(3,2);
t661 = sin(t685);
t664 = cos(t685);
t551 = t602 * t664 + (t581 * t655 + t584 * t649) * t661;
t554 = -t649 * t581 + t584 * t655;
t751 = t678 * t689;
t614 = t675 * t751 - t677 * t695;
t617 = t675 * t695 + t677 * t751;
t712 = t614 * t655 + t649 * t617;
t563 = t712 * t661 + t664 * t757;
t578 = -t649 * t614 + t617 * t655;
t611 = t677 * t649 + t675 * t655;
t587 = -t611 * t676 * t661 + t678 * t664;
t682 = legFrame(3,1);
t652 = sin(t682);
t658 = cos(t682);
t608 = -t675 * t649 + t655 * t677;
t776 = t608 * t658;
t530 = -(t563 * t652 - t578 * t658) * t823 + (-t551 * t652 + t554 * t658) * t694 - (t587 * t652 + t676 * t776) * t819;
t575 = 0.1e1 / (pkin(2) * t748 + t602 * t694 + t757 * t823);
t812 = t530 * t575;
t769 = t652 * t608;
t531 = (t563 * t658 + t578 * t652) * t823 + (t551 * t658 + t652 * t554) * t694 + (t587 * t658 - t676 * t769) * t819;
t811 = t531 * t575;
t697 = cos(qJ(2,2));
t745 = t691 * t702;
t633 = pkin(2) * t697 + t745;
t646 = t702 * t697;
t630 = pkin(2) * t691 - t646;
t708 = -t630 * t678 + t690 * t820;
t582 = -t633 * t677 - t708 * t675;
t585 = t633 * t675 - t708 * t677;
t746 = t690 * t678;
t603 = pkin(3) * t746 + t630 * t676;
t680 = legFrame(2,3);
t650 = sin(t680);
t656 = cos(t680);
t686 = legFrame(2,2);
t662 = sin(t686);
t665 = cos(t686);
t552 = t603 * t665 + (t582 * t656 + t585 * t650) * t662;
t555 = -t650 * t582 + t585 * t656;
t750 = t678 * t691;
t615 = t675 * t750 - t677 * t697;
t618 = t675 * t697 + t677 * t750;
t711 = t615 * t656 + t650 * t618;
t564 = t711 * t662 + t665 * t756;
t579 = -t650 * t615 + t618 * t656;
t612 = t677 * t650 + t675 * t656;
t588 = -t612 * t676 * t662 + t678 * t665;
t683 = legFrame(2,1);
t653 = sin(t683);
t659 = cos(t683);
t609 = -t675 * t650 + t656 * t677;
t775 = t609 * t659;
t532 = -(t564 * t653 - t579 * t659) * t822 + (-t552 * t653 + t555 * t659) * t696 - (t588 * t653 + t676 * t775) * t818;
t576 = 0.1e1 / (pkin(2) * t746 + t603 * t696 + t756 * t822);
t810 = t532 * t576;
t767 = t653 * t609;
t533 = (t564 * t659 + t579 * t653) * t822 + (t552 * t659 + t653 * t555) * t696 + (t588 * t659 - t676 * t767) * t818;
t809 = t533 * t576;
t699 = cos(qJ(2,1));
t743 = t693 * t702;
t634 = pkin(2) * t699 + t743;
t647 = t702 * t699;
t631 = pkin(2) * t693 - t647;
t707 = -t631 * t678 + t692 * t820;
t583 = -t634 * t677 - t707 * t675;
t586 = t634 * t675 - t707 * t677;
t744 = t692 * t678;
t604 = pkin(3) * t744 + t631 * t676;
t681 = legFrame(1,3);
t651 = sin(t681);
t657 = cos(t681);
t687 = legFrame(1,2);
t663 = sin(t687);
t666 = cos(t687);
t553 = t604 * t666 + (t583 * t657 + t586 * t651) * t663;
t556 = -t651 * t583 + t586 * t657;
t749 = t678 * t693;
t616 = t675 * t749 - t677 * t699;
t619 = t675 * t699 + t677 * t749;
t710 = t616 * t657 + t651 * t619;
t565 = t710 * t663 + t666 * t755;
t580 = -t651 * t616 + t619 * t657;
t613 = t677 * t651 + t675 * t657;
t589 = -t613 * t676 * t663 + t678 * t666;
t684 = legFrame(1,1);
t654 = sin(t684);
t660 = cos(t684);
t610 = -t675 * t651 + t657 * t677;
t774 = t610 * t660;
t534 = -(t565 * t654 - t580 * t660) * t821 + (-t553 * t654 + t556 * t660) * t698 - (t589 * t654 + t676 * t774) * t817;
t577 = 0.1e1 / (pkin(2) * t744 + t604 * t698 + t755 * t821);
t808 = t534 * t577;
t765 = t654 * t610;
t535 = (t565 * t660 + t580 * t654) * t821 + (t553 * t660 + t654 * t556) * t698 + (t589 * t660 - t676 * t765) * t817;
t807 = t535 * t577;
t763 = t658 * t661;
t569 = t608 * t763 - t652 * t611;
t754 = t676 * t694;
t539 = (t578 * t763 - t652 * t712) * t688 + t569 * t754;
t641 = pkin(2) + t816;
t626 = t641 * t748;
t590 = 0.1e1 / (t626 + (t689 * t816 + t629) * t754);
t806 = t539 * t590;
t762 = t659 * t662;
t570 = t609 * t762 - t653 * t612;
t753 = t676 * t696;
t540 = (t579 * t762 - t653 * t711) * t690 + t570 * t753;
t642 = pkin(2) + t815;
t627 = t642 * t746;
t591 = 0.1e1 / (t627 + (t691 * t815 + t630) * t753);
t805 = t540 * t591;
t761 = t660 * t663;
t571 = t610 * t761 - t654 * t613;
t752 = t676 * t698;
t541 = (t580 * t761 - t654 * t710) * t692 + t571 * t752;
t643 = pkin(2) + t814;
t628 = t643 * t744;
t592 = 0.1e1 / (t628 + (t693 * t814 + t631) * t752);
t804 = t541 * t592;
t768 = t652 * t661;
t572 = t608 * t768 + t611 * t658;
t542 = (-t578 * t768 - t712 * t658) * t688 - t572 * t754;
t803 = t542 * t590;
t766 = t653 * t662;
t573 = t609 * t766 + t612 * t659;
t543 = (-t579 * t766 - t711 * t659) * t690 - t573 * t753;
t802 = t543 * t591;
t764 = t654 * t663;
t574 = t610 * t764 + t613 * t660;
t544 = (-t580 * t764 - t710 * t660) * t692 - t574 * t752;
t801 = t544 * t592;
t623 = t689 * t641 - t645;
t773 = (t641 * t695 + t747) * t678;
t545 = -t572 * t773 + (t611 * t768 - t776) * t623;
t593 = 0.1e1 / (t623 * t754 + t626);
t800 = t545 * t593;
t624 = t691 * t642 - t646;
t772 = (t642 * t697 + t745) * t678;
t546 = -t573 * t772 + (t612 * t766 - t775) * t624;
t594 = 0.1e1 / (t624 * t753 + t627);
t799 = t546 * t594;
t625 = t693 * t643 - t647;
t771 = (t643 * t699 + t743) * t678;
t547 = -t574 * t771 + (t613 * t764 - t774) * t625;
t595 = 0.1e1 / (t625 * t752 + t628);
t798 = t547 * t595;
t548 = t569 * t773 - (t611 * t763 + t769) * t623;
t797 = t548 * t593;
t549 = t570 * t772 - (t612 * t762 + t767) * t624;
t796 = t549 * t594;
t550 = t571 * t771 - (t613 * t761 + t765) * t625;
t795 = t550 * t595;
t557 = t608 * t754 + t688 * (t608 * t751 + t695 * t611);
t794 = t557 * t664;
t558 = t609 * t753 + t690 * (t609 * t750 + t697 * t612);
t793 = t558 * t665;
t559 = t610 * t752 + t692 * (t610 * t749 + t699 * t613);
t792 = t559 * t666;
t704 = rSges(3,2) ^ 2;
t705 = rSges(3,1) ^ 2;
t635 = (-t704 + t705) * m(3) + Icges(3,2) - Icges(3,1);
t667 = 2 * rSges(3,1) * t703;
t716 = Icges(3,1) + Icges(2,3) + (pkin(2) ^ 2 + t700 ^ 2 + t704) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t738 = -0.2e1 * pkin(2) * t827;
t566 = t635 * t672 + (t688 * t828 + t667) * t694 + t688 * t738 + t716;
t791 = t566 * t590;
t567 = t635 * t673 + (t690 * t828 + t667) * t696 + t690 * t738 + t716;
t790 = t567 * t591;
t568 = t635 * t674 + (t692 * t828 + t667) * t698 + t692 * t738 + t716;
t789 = t568 * t592;
t671 = m(1) + m(2) + m(3);
t788 = t575 * t671;
t787 = t576 * t671;
t786 = t577 * t671;
t639 = -rSges(3,2) * t813 + Icges(3,6);
t640 = rSges(3,1) * t813 - Icges(3,5);
t605 = t639 * t694 - t688 * t640;
t785 = t590 * t605;
t606 = t639 * t696 - t690 * t640;
t784 = t591 * t606;
t607 = t639 * t698 - t692 * t640;
t783 = t592 * t607;
t706 = 0.1e1 / pkin(3);
t782 = t593 * t706;
t781 = t594 * t706;
t780 = t595 * t706;
t638 = m(2) * rSges(2,2) - t813;
t742 = m(2) * rSges(2,1) + t703;
t596 = (t715 * m(3) + t742) * t695 - t689 * t638;
t779 = t596 * t676;
t597 = (t714 * m(3) + t742) * t697 - t691 * t638;
t778 = t597 * t676;
t598 = (t713 * m(3) + t742) * t699 - t693 * t638;
t777 = t598 * t676;
t637 = (t704 + t705) * m(3) + Icges(3,3);
t770 = t637 * t706;
t760 = t664 * t676;
t759 = t665 * t676;
t758 = t666 * t676;
t741 = t575 * t826;
t740 = t576 * t825;
t739 = t577 * t824;
t737 = (-t608 * t773 + t611 * t623) * t593 * t664;
t736 = (-t609 * t772 + t612 * t624) * t594 * t665;
t735 = (-t610 * t771 + t613 * t625) * t595 * t666;
t734 = t575 * t779;
t733 = t576 * t778;
t732 = t577 * t777;
t731 = t590 * t779;
t730 = t591 * t778;
t729 = t592 * t777;
t728 = t605 * t782;
t727 = t593 * t770;
t726 = t606 * t781;
t725 = t594 * t770;
t724 = t607 * t780;
t723 = t595 * t770;
t722 = t782 * t826;
t721 = t781 * t825;
t720 = t780 * t824;
t719 = t706 * t737;
t718 = t706 * t736;
t717 = t706 * t735;
t538 = -((-t699 * t610 + t613 * t749) * t666 - t663 * t755) * t821 + ((t634 * t610 + t707 * t613) * t666 + t604 * t663) * t698 + (t613 * t758 + t678 * t663) * t817;
t537 = -((-t697 * t609 + t612 * t750) * t665 - t662 * t756) * t822 + ((t633 * t609 + t708 * t612) * t665 + t603 * t662) * t696 + (t612 * t759 + t678 * t662) * t818;
t536 = -((-t695 * t608 + t611 * t751) * t664 - t661 * t757) * t823 + ((t632 * t608 + t709 * t611) * t664 + t602 * t661) * t694 + (t611 * t760 + t678 * t661) * t819;
t529 = t637 * t717 + (t538 * t824 - t607 * t792) * t577;
t528 = t637 * t718 + (t537 * t825 - t606 * t793) * t576;
t527 = t637 * t719 + (t536 * t826 - t605 * t794) * t575;
t526 = t717 * t824 + (-t559 * t598 * t758 + t538 * t671) * t577;
t525 = t718 * t825 + (-t558 * t597 * t759 + t537 * t671) * t576;
t524 = t719 * t826 + (-t557 * t596 * t760 + t536 * t671) * t575;
t523 = t607 * t717 + (t538 * t777 - t568 * t792) * t577;
t522 = t606 * t718 + (t537 * t778 - t567 * t793) * t576;
t521 = t605 * t719 + (t536 * t779 - t566 * t794) * t575;
t520 = t535 * t739 + t541 * t783 + t550 * t723;
t519 = t534 * t739 + t544 * t783 + t547 * t723;
t518 = t533 * t740 + t540 * t784 + t549 * t725;
t517 = t532 * t740 + t543 * t784 + t546 * t725;
t516 = t531 * t741 + t539 * t785 + t548 * t727;
t515 = t530 * t741 + t542 * t785 + t545 * t727;
t514 = t535 * t786 + t541 * t729 + t550 * t720;
t513 = t534 * t786 + t544 * t729 + t547 * t720;
t512 = t533 * t787 + t540 * t730 + t549 * t721;
t511 = t532 * t787 + t543 * t730 + t546 * t721;
t510 = t531 * t788 + t539 * t731 + t548 * t722;
t509 = t530 * t788 + t542 * t731 + t545 * t722;
t508 = t535 * t732 + t541 * t789 + t550 * t724;
t507 = t534 * t732 + t544 * t789 + t547 * t724;
t506 = t533 * t733 + t540 * t790 + t549 * t726;
t505 = t532 * t733 + t543 * t790 + t546 * t726;
t504 = t531 * t734 + t539 * t791 + t548 * t728;
t503 = t530 * t734 + t542 * t791 + t545 * t728;
t1 = [m(4) + (-t523 * t792 + t526 * t538) * t577 + (-t522 * t793 + t525 * t537) * t576 + (-t521 * t794 + t524 * t536) * t575 + (t527 * t737 + t528 * t736 + t529 * t735) * t706, t521 * t803 + t522 * t802 + t523 * t801 + t524 * t812 + t525 * t810 + t526 * t808 + (t527 * t800 + t528 * t799 + t529 * t798) * t706, t521 * t806 + t522 * t805 + t523 * t804 + t524 * t811 + t525 * t809 + t526 * t807 + (t527 * t797 + t528 * t796 + t529 * t795) * t706; (-t507 * t792 + t513 * t538) * t577 + (-t505 * t793 + t511 * t537) * t576 + (-t503 * t794 + t509 * t536) * t575 + (t515 * t737 + t517 * t736 + t519 * t735) * t706, t503 * t803 + t505 * t802 + t507 * t801 + t509 * t812 + t511 * t810 + t513 * t808 + m(4) + (t515 * t800 + t517 * t799 + t519 * t798) * t706, t503 * t806 + t505 * t805 + t507 * t804 + t509 * t811 + t511 * t809 + t513 * t807 + (t515 * t797 + t517 * t796 + t519 * t795) * t706; (-t508 * t792 + t514 * t538) * t577 + (-t506 * t793 + t512 * t537) * t576 + (-t504 * t794 + t510 * t536) * t575 + (t516 * t737 + t518 * t736 + t520 * t735) * t706, t504 * t803 + t506 * t802 + t508 * t801 + t510 * t812 + t512 * t810 + t514 * t808 + (t516 * t800 + t518 * t799 + t520 * t798) * t706, t504 * t806 + t506 * t805 + t508 * t804 + t510 * t811 + t512 * t809 + t514 * t807 + m(4) + (t516 * t797 + t518 * t796 + t520 * t795) * t706;];
MX  = t1;
