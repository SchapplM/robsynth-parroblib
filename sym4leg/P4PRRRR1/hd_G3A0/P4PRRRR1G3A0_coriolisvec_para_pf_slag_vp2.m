% Calculate vector of centrifugal and coriolis load on the joints for
% P4PRRRR1G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [4x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:02:59
% EndTime: 2020-03-02 19:03:07
% DurationCPUTime: 7.80s
% Computational Cost: add. (3340->324), mult. (6504->702), div. (3852->18), fcn. (6932->26), ass. (0->280)
t682 = cos(qJ(3,4));
t658 = t682 ^ 2;
t659 = 0.1e1 / t682;
t661 = t659 / t658;
t680 = sin(qJ(3,4));
t840 = t661 * t680;
t696 = cos(qJ(3,3));
t666 = t696 ^ 2;
t667 = 0.1e1 / t696;
t669 = t667 / t666;
t690 = sin(qJ(3,3));
t839 = t669 * t690;
t698 = cos(qJ(3,2));
t670 = t698 ^ 2;
t671 = 0.1e1 / t698;
t673 = t671 / t670;
t692 = sin(qJ(3,2));
t838 = t673 * t692;
t700 = cos(qJ(3,1));
t674 = t700 ^ 2;
t675 = 0.1e1 / t700;
t677 = t675 / t674;
t694 = sin(qJ(3,1));
t837 = t677 * t694;
t706 = xP(4);
t655 = sin(t706);
t656 = cos(t706);
t709 = koppelP(4,2);
t713 = koppelP(4,1);
t626 = t655 * t713 + t656 * t709;
t630 = -t655 * t709 + t656 * t713;
t686 = legFrame(4,2);
t647 = sin(t686);
t651 = cos(t686);
t681 = sin(qJ(2,4));
t657 = 0.1e1 / t681;
t660 = 0.1e1 / t682 ^ 2;
t683 = cos(qJ(2,4));
t702 = xDP(4);
t704 = xDP(2);
t705 = xDP(1);
t717 = 0.1e1 / pkin(2);
t703 = xDP(3);
t777 = t680 * t703;
t571 = (-t683 * t660 * t777 + (-t651 * (-t626 * t702 + t705) + t647 * (t630 * t702 + t704)) * t659) * t717 * t657;
t570 = t571 ^ 2;
t770 = t703 * t717;
t646 = mrSges(3,2) * t770;
t730 = -mrSges(3,1) * t682 + mrSges(3,2) * t680;
t758 = t659 * t770;
t679 = t703 ^ 2;
t718 = 0.1e1 / pkin(2) ^ 2;
t780 = t679 * t718;
t817 = t571 * t683;
t821 = mrSges(3,1) * t680;
t685 = mrSges(2,2) - mrSges(3,3);
t823 = t685 / 0.2e1;
t534 = (-mrSges(2,1) * t570 + t730 * (t660 * t780 + t570)) * t681 - 0.2e1 * (t571 * t823 + t758 * t821 + t646) * t817;
t779 = t680 * t682;
t803 = t660 * t717;
t546 = ((pkin(2) * t658 * t817 - t681 * t777) * t571 * t803 + (-t571 * t681 * t779 + t683 * t758) * t661 * t770) * t657;
t781 = t679 * t717;
t562 = (-pkin(2) * t570 * t682 - t661 * t781) * t657;
t610 = (mrSges(2,1) - t730) * t683 - t681 * t685;
t662 = m(1) + m(2) + m(3);
t740 = t780 * t840;
t813 = (mrSges(3,2) * t682 + t821) * t681;
t539 = -t610 * t546 - t662 * t562 - t740 * t813;
t836 = t539 + t534;
t710 = koppelP(3,2);
t714 = koppelP(3,1);
t627 = t655 * t714 + t656 * t710;
t631 = -t655 * t710 + t656 * t714;
t687 = legFrame(3,2);
t648 = sin(t687);
t652 = cos(t687);
t691 = sin(qJ(2,3));
t663 = 0.1e1 / t691;
t668 = 0.1e1 / t696 ^ 2;
t697 = cos(qJ(2,3));
t775 = t690 * t703;
t575 = (-t697 * t668 * t775 + (-t652 * (-t627 * t702 + t705) + t648 * (t631 * t702 + t704)) * t667) * t717 * t663;
t572 = t575 ^ 2;
t729 = -mrSges(3,1) * t696 + mrSges(3,2) * t690;
t750 = t667 * t770;
t816 = t575 * t697;
t820 = mrSges(3,1) * t690;
t536 = (-mrSges(2,1) * t572 + t729 * (t668 * t780 + t572)) * t691 - 0.2e1 * (t575 * t823 + t750 * t820 + t646) * t816;
t776 = t690 * t696;
t788 = t668 * t717;
t547 = ((pkin(2) * t666 * t816 - t691 * t775) * t575 * t788 + (-t575 * t691 * t776 + t697 * t750) * t669 * t770) * t663;
t563 = (-pkin(2) * t572 * t696 - t669 * t781) * t663;
t611 = (mrSges(2,1) - t729) * t697 - t691 * t685;
t739 = t780 * t839;
t812 = (mrSges(3,2) * t696 + t820) * t691;
t543 = -t611 * t547 - t662 * t563 - t739 * t812;
t835 = t543 + t536;
t711 = koppelP(2,2);
t715 = koppelP(2,1);
t628 = t655 * t715 + t656 * t711;
t632 = -t655 * t711 + t656 * t715;
t688 = legFrame(2,2);
t649 = sin(t688);
t653 = cos(t688);
t693 = sin(qJ(2,2));
t664 = 0.1e1 / t693;
t672 = 0.1e1 / t698 ^ 2;
t699 = cos(qJ(2,2));
t773 = t692 * t703;
t576 = (-t699 * t672 * t773 + (-t653 * (-t628 * t702 + t705) + t649 * (t632 * t702 + t704)) * t671) * t717 * t664;
t573 = t576 ^ 2;
t728 = -mrSges(3,1) * t698 + mrSges(3,2) * t692;
t748 = t671 * t770;
t815 = t576 * t699;
t819 = mrSges(3,1) * t692;
t537 = (-mrSges(2,1) * t573 + t728 * (t672 * t780 + t573)) * t693 - 0.2e1 * (t576 * t823 + t748 * t819 + t646) * t815;
t774 = t692 * t698;
t785 = t672 * t717;
t548 = ((pkin(2) * t670 * t815 - t693 * t773) * t576 * t785 + (-t576 * t693 * t774 + t699 * t748) * t673 * t770) * t664;
t564 = (-pkin(2) * t573 * t698 - t673 * t781) * t664;
t612 = (mrSges(2,1) - t728) * t699 - t693 * t685;
t738 = t780 * t838;
t811 = (mrSges(3,2) * t698 + t819) * t693;
t544 = -t612 * t548 - t662 * t564 - t738 * t811;
t834 = t544 + t537;
t712 = koppelP(1,2);
t716 = koppelP(1,1);
t629 = t655 * t716 + t656 * t712;
t633 = -t655 * t712 + t656 * t716;
t689 = legFrame(1,2);
t650 = sin(t689);
t654 = cos(t689);
t695 = sin(qJ(2,1));
t665 = 0.1e1 / t695;
t676 = 0.1e1 / t700 ^ 2;
t701 = cos(qJ(2,1));
t771 = t694 * t703;
t577 = (-t701 * t676 * t771 + (-t654 * (-t629 * t702 + t705) + t650 * (t633 * t702 + t704)) * t675) * t717 * t665;
t574 = t577 ^ 2;
t727 = -mrSges(3,1) * t700 + mrSges(3,2) * t694;
t746 = t675 * t770;
t814 = t577 * t701;
t818 = mrSges(3,1) * t694;
t538 = (-mrSges(2,1) * t574 + t727 * (t676 * t780 + t574)) * t695 - 0.2e1 * (t577 * t823 + t746 * t818 + t646) * t814;
t772 = t694 * t700;
t782 = t676 * t717;
t549 = ((pkin(2) * t674 * t814 - t695 * t771) * t577 * t782 + (-t577 * t695 * t772 + t701 * t746) * t677 * t770) * t665;
t565 = (-pkin(2) * t574 * t700 - t677 * t781) * t665;
t613 = (mrSges(2,1) - t727) * t701 - t695 * t685;
t737 = t780 * t837;
t810 = (mrSges(3,2) * t700 + t818) * t695;
t545 = -t613 * t549 - t662 * t565 - t737 * t810;
t833 = t545 + t538;
t684 = Ifges(3,1) - Ifges(3,2);
t822 = Ifges(3,1) + Ifges(2,3);
t828 = 2 * Ifges(3,4);
t619 = -t674 * t684 + t772 * t828 + t822;
t638 = Ifges(3,5) * t694 + Ifges(3,6) * t700;
t542 = -t619 * t549 - t613 * t565 + t638 * t737;
t731 = -Ifges(3,6) * t770 / 0.2e1;
t732 = Ifges(3,5) * t770 / 0.2e1;
t784 = t675 * t694;
t824 = 0.2e1 * t674;
t553 = (t577 * t684 * t694 + t675 * t732) * t700 + t731 * t784 + (t824 - 0.1e1) * t577 * Ifges(3,4);
t769 = t703 * t718;
t741 = t553 * t665 * t769;
t783 = t675 * t717;
t751 = t665 * t783;
t832 = t542 * t751 + 0.2e1 * t676 * t741;
t618 = -t670 * t684 + t774 * t828 + t822;
t637 = Ifges(3,5) * t692 + Ifges(3,6) * t698;
t541 = -t618 * t548 - t612 * t564 + t637 * t738;
t787 = t671 * t692;
t825 = 0.2e1 * t670;
t552 = (t576 * t684 * t692 + t671 * t732) * t698 + t731 * t787 + (t825 - 0.1e1) * t576 * Ifges(3,4);
t742 = t552 * t664 * t769;
t786 = t671 * t717;
t753 = t664 * t786;
t831 = t541 * t753 + 0.2e1 * t672 * t742;
t617 = -t666 * t684 + t776 * t828 + t822;
t636 = Ifges(3,5) * t690 + Ifges(3,6) * t696;
t540 = -t617 * t547 - t611 * t563 + t636 * t739;
t790 = t667 * t690;
t826 = 0.2e1 * t666;
t551 = (t575 * t684 * t690 + t667 * t732) * t696 + t731 * t790 + (t826 - 0.1e1) * t575 * Ifges(3,4);
t743 = t551 * t663 * t769;
t789 = t667 * t717;
t755 = t663 * t789;
t830 = t540 * t755 + 0.2e1 * t668 * t743;
t614 = -t658 * t684 + t779 * t828 + t822;
t634 = Ifges(3,5) * t680 + Ifges(3,6) * t682;
t535 = -t614 * t546 - t610 * t562 + t634 * t740;
t778 = t680 * t684;
t805 = t659 * t680;
t827 = 0.2e1 * t658;
t550 = (t571 * t778 + t659 * t732) * t682 + t731 * t805 + (t827 - 0.1e1) * t571 * Ifges(3,4);
t744 = t550 * t657 * t769;
t804 = t659 * t717;
t759 = t657 * t804;
t829 = t535 * t759 + 0.2e1 * t660 * t744;
t678 = t702 ^ 2;
t809 = t657 * t534;
t808 = t657 * t539;
t807 = t657 * t680;
t806 = t657 * t678;
t802 = t663 * t536;
t801 = t663 * t543;
t800 = t663 * t690;
t799 = t663 * t678;
t798 = t664 * t537;
t797 = t664 * t544;
t796 = t664 * t692;
t795 = t664 * t678;
t794 = t665 * t538;
t793 = t665 * t545;
t792 = t665 * t694;
t791 = t665 * t678;
t768 = t647 * t804;
t767 = t648 * t789;
t766 = t649 * t786;
t765 = t650 * t783;
t764 = t651 * t804;
t763 = t652 * t789;
t762 = t653 * t786;
t761 = t654 * t783;
t757 = t683 * t803;
t749 = t697 * t788;
t747 = t699 * t785;
t745 = t701 * t782;
t708 = mrSges(4,1);
t707 = mrSges(4,2);
t625 = t650 * t695 + t654 * t701;
t624 = -t650 * t701 + t654 * t695;
t623 = t649 * t693 + t653 * t699;
t622 = -t649 * t699 + t653 * t693;
t621 = t648 * t691 + t652 * t697;
t620 = -t648 * t697 + t652 * t691;
t616 = t647 * t681 + t651 * t683;
t615 = -t647 * t683 + t651 * t681;
t609 = (t629 * t654 + t633 * t650) * t751;
t608 = (t628 * t653 + t632 * t649) * t753;
t607 = (t627 * t652 + t631 * t648) * t755;
t606 = (t626 * t651 + t630 * t647) * t759;
t605 = (-t613 * t761 + t625 * t662) * t665;
t604 = (t613 * t765 + t624 * t662) * t665;
t603 = (-t612 * t762 + t623 * t662) * t664;
t602 = (t612 * t766 + t622 * t662) * t664;
t601 = (-t611 * t763 + t621 * t662) * t663;
t600 = (t611 * t767 + t620 * t662) * t663;
t599 = (-t610 * t764 + t616 * t662) * t657;
t598 = (t610 * t768 + t615 * t662) * t657;
t597 = -t783 * t810 + (-t613 * t745 + t662 * t675) * t792;
t596 = -t786 * t811 + (-t612 * t747 + t662 * t671) * t796;
t595 = -t789 * t812 + (-t611 * t749 + t662 * t667) * t800;
t594 = -t804 * t813 + (-t610 * t757 + t659 * t662) * t807;
t593 = (t624 * t633 - t625 * t629) * t665;
t592 = (t622 * t632 - t623 * t628) * t664;
t591 = (t620 * t631 - t621 * t627) * t663;
t590 = (t615 * t630 - t616 * t626) * t657;
t589 = (t613 * t625 - t619 * t761) * t665;
t588 = (t613 * t624 + t619 * t765) * t665;
t587 = (t612 * t623 - t618 * t762) * t664;
t586 = (t612 * t622 + t618 * t766) * t664;
t585 = (t611 * t621 - t617 * t763) * t663;
t584 = (t611 * t620 + t617 * t767) * t663;
t583 = (t610 * t616 - t614 * t764) * t657;
t582 = (t610 * t615 + t614 * t768) * t657;
t581 = t638 * t783 + (t613 * t675 - t619 * t745) * t792;
t580 = t637 * t786 + (t612 * t671 - t618 * t747) * t796;
t579 = t636 * t789 + (t611 * t667 - t617 * t749) * t800;
t578 = t634 * t804 + (t610 * t659 - t614 * t757) * t807;
t561 = t593 * t662 + t609 * t613;
t560 = t592 * t662 + t608 * t612;
t559 = t591 * t662 + t607 * t611;
t558 = t590 * t662 + t606 * t610;
t557 = t593 * t613 + t609 * t619;
t556 = t592 * t612 + t608 * t618;
t555 = t591 * t611 + t607 * t617;
t554 = t590 * t610 + t606 * t614;
t1 = [(-(-t589 * t761 + t605 * t625) * t633 - (t589 * t765 + t605 * t624) * t629) * t791 + t625 * t793 + t625 * t794 + (-(-t587 * t762 + t603 * t623) * t632 - (t587 * t766 + t603 * t622) * t628) * t795 + t623 * t797 + t623 * t798 + (-(-t585 * t763 + t601 * t621) * t631 - (t585 * t767 + t601 * t620) * t627) * t799 + t621 * t801 + t621 * t802 + (-(-t583 * t764 + t599 * t616) * t630 - (t583 * t768 + t599 * t615) * t626) * t806 + t616 * t808 + t616 * t809 + t678 * (t655 * t707 - t656 * t708) - t832 * t654 - t831 * t653 - t830 * t652 - t829 * t651; (-(-t588 * t761 + t604 * t625) * t633 - (t588 * t765 + t604 * t624) * t629) * t791 + t624 * t793 + t624 * t794 + (-(-t586 * t762 + t602 * t623) * t632 - (t586 * t766 + t602 * t622) * t628) * t795 + t622 * t797 + t622 * t798 + (-(-t584 * t763 + t600 * t621) * t631 - (t584 * t767 + t600 * t620) * t627) * t799 + t620 * t801 + t620 * t802 + (-(-t582 * t764 + t598 * t616) * t630 - (t582 * t768 + t598 * t615) * t626) * t806 + t615 * t808 + t615 * t809 - t678 * (t655 * t708 + t656 * t707) + t832 * t650 + t831 * t649 + t830 * t648 + t829 * t647; (-(-t579 * t763 + t595 * t621) * t631 - (t579 * t767 + t595 * t620) * t627) * t799 + (-(-t578 * t764 + t594 * t616) * t630 - (t578 * t768 + t594 * t615) * t626) * t806 + (-(-t581 * t761 + t597 * t625) * t633 - (t581 * t765 + t597 * t624) * t629) * t791 + (-(-t580 * t762 + t596 * t623) * t632 - (t580 * t766 + t596 * t622) * t628) * t795 - t535 * t757 * t807 - t542 * t745 * t792 - t541 * t747 * t796 - t540 * t749 * t800 - 0.2e1 * t697 * t743 * t839 - 0.2e1 * t683 * t744 * t840 - 0.2e1 * t701 * t741 * t837 - 0.2e1 * t699 * t742 * t838 + (-t570 * (Ifges(3,4) * t827 + t682 * t778 - Ifges(3,4)) + Ifges(3,3) * t740 - t634 * t546 + t562 * t813) * t804 + (-t572 * (Ifges(3,4) * t826 + t684 * t776 - Ifges(3,4)) + Ifges(3,3) * t739 - t636 * t547 + t563 * t812) * t789 + (-t573 * (Ifges(3,4) * t825 + t684 * t774 - Ifges(3,4)) + Ifges(3,3) * t738 - t637 * t548 + t564 * t811) * t786 + (-t574 * (Ifges(3,4) * t824 + t684 * t772 - Ifges(3,4)) + Ifges(3,3) * t737 - t638 * t549 + t565 * t810) * t783 + t836 * t657 * t805 + t835 * t663 * t790 + t834 * t664 * t787 + t833 * t665 * t784; (-(-t557 * t761 + t561 * t625) * t633 - (t557 * t765 + t561 * t624) * t629) * t791 + (-(-t556 * t762 + t560 * t623) * t632 - (t556 * t766 + t560 * t622) * t628) * t795 + (-(-t555 * t763 + t559 * t621) * t631 - (t555 * t767 + t559 * t620) * t627) * t799 + (-(-t554 * t764 + t558 * t616) * t630 - (t554 * t768 + t558 * t615) * t626) * t806 + (0.2e1 * t553 * t746 + t542) * t609 + (0.2e1 * t552 * t748 + t541) * t608 + (0.2e1 * t551 * t750 + t540) * t607 + (0.2e1 * t550 * t758 + t535) * t606 + t833 * t593 + t834 * t592 + t835 * t591 + t836 * t590;];
taucX  = t1;
