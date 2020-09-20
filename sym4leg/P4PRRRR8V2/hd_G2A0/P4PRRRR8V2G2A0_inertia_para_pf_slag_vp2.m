% Calculate inertia matrix for parallel robot
% P4PRRRR8V2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:22
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:16:56
% EndTime: 2020-08-07 11:17:00
% DurationCPUTime: 3.97s
% Computational Cost: add. (15484->458), mult. (32520->856), div. (1280->5), fcn. (30920->30), ass. (0->304)
t735 = sin(qJ(3,4));
t882 = pkin(2) * t735;
t744 = sin(qJ(3,3));
t881 = pkin(2) * t744;
t746 = sin(qJ(3,2));
t880 = pkin(2) * t746;
t748 = sin(qJ(3,1));
t879 = pkin(2) * t748;
t737 = cos(qJ(3,4));
t725 = t737 ^ 2;
t878 = pkin(3) * t725;
t750 = cos(qJ(3,3));
t728 = t750 ^ 2;
t877 = pkin(3) * t728;
t752 = cos(qJ(3,2));
t729 = t752 ^ 2;
t876 = pkin(3) * t729;
t754 = cos(qJ(3,1));
t730 = t754 ^ 2;
t875 = pkin(3) * t730;
t874 = pkin(3) * t737;
t873 = pkin(3) * t750;
t872 = pkin(3) * t752;
t871 = pkin(3) * t754;
t772 = 0.1e1 / pkin(3);
t870 = Ifges(3,3) * t772;
t869 = m(3) * pkin(2) + mrSges(2,1);
t736 = sin(qJ(2,4));
t738 = cos(qJ(2,4));
t757 = pkin(7) + pkin(6);
t699 = pkin(2) * t736 - t757 * t738;
t700 = pkin(2) * t738 + t736 * t757;
t731 = sin(pkin(8));
t733 = cos(pkin(8));
t734 = cos(pkin(4));
t819 = t734 * t738;
t825 = t733 * t734;
t639 = (t731 * t736 - t733 * t819) * t874 - t700 * t825 + t699 * t731;
t732 = sin(pkin(4));
t812 = t735 * t734;
t675 = pkin(3) * t812 + t699 * t732;
t832 = t732 * t736;
t647 = 0.1e1 / (pkin(2) * t812 + t675 * t737 + t832 * t878);
t868 = t639 * t647;
t835 = t731 * t734;
t640 = (t731 * t819 + t733 * t736) * t874 + t700 * t835 + t699 * t733;
t867 = t640 * t647;
t866 = t640 * t772;
t745 = sin(qJ(2,3));
t751 = cos(qJ(2,3));
t701 = pkin(2) * t745 - t757 * t751;
t704 = pkin(2) * t751 + t745 * t757;
t815 = t734 * t751;
t641 = (t731 * t745 - t733 * t815) * t873 - t704 * t825 + t701 * t731;
t811 = t744 * t734;
t676 = pkin(3) * t811 + t701 * t732;
t830 = t732 * t745;
t648 = 0.1e1 / (pkin(2) * t811 + t676 * t750 + t830 * t877);
t865 = t641 * t648;
t747 = sin(qJ(2,2));
t753 = cos(qJ(2,2));
t702 = pkin(2) * t747 - t757 * t753;
t705 = pkin(2) * t753 + t747 * t757;
t814 = t734 * t753;
t642 = (t731 * t747 - t733 * t814) * t872 - t705 * t825 + t702 * t731;
t810 = t746 * t734;
t677 = pkin(3) * t810 + t702 * t732;
t828 = t732 * t747;
t649 = 0.1e1 / (pkin(2) * t810 + t677 * t752 + t828 * t876);
t864 = t642 * t649;
t749 = sin(qJ(2,1));
t755 = cos(qJ(2,1));
t703 = pkin(2) * t749 - t757 * t755;
t706 = pkin(2) * t755 + t749 * t757;
t813 = t734 * t755;
t643 = (t731 * t749 - t733 * t813) * t871 - t706 * t825 + t703 * t731;
t809 = t748 * t734;
t678 = pkin(3) * t809 + t703 * t732;
t826 = t732 * t749;
t650 = 0.1e1 / (pkin(2) * t809 + t678 * t754 + t826 * t875);
t863 = t643 * t650;
t644 = (t731 * t815 + t733 * t745) * t873 + t704 * t835 + t701 * t733;
t862 = t644 * t648;
t861 = t644 * t772;
t645 = (t731 * t814 + t733 * t747) * t872 + t705 * t835 + t702 * t733;
t860 = t645 * t649;
t859 = t645 * t772;
t646 = (t731 * t813 + t733 * t749) * t871 + t706 * t835 + t703 * t733;
t858 = t646 * t650;
t857 = t646 * t772;
t799 = t737 * mrSges(3,1) - mrSges(3,2) * t735;
t659 = t799 * t734 - (t735 * mrSges(3,1) + t737 * mrSges(3,2)) * t832;
t856 = t659 * t772;
t820 = t734 * t736;
t679 = t731 * t820 - t733 * t738;
t836 = t731 * t732;
t660 = t735 * t679 + t737 * t836;
t740 = legFrame(4,2);
t715 = sin(t740);
t855 = t660 * t715;
t719 = cos(t740);
t854 = t660 * t719;
t818 = t734 * t745;
t682 = t731 * t818 - t733 * t751;
t662 = t744 * t682 + t750 * t836;
t741 = legFrame(3,2);
t716 = sin(t741);
t853 = t662 * t716;
t720 = cos(t741);
t852 = t662 * t720;
t817 = t734 * t747;
t683 = t731 * t817 - t733 * t753;
t663 = t746 * t683 + t752 * t836;
t742 = legFrame(2,2);
t717 = sin(t742);
t851 = t663 * t717;
t721 = cos(t742);
t850 = t663 * t721;
t816 = t734 * t749;
t684 = t731 * t816 - t733 * t755;
t664 = t748 * t684 + t754 * t836;
t743 = legFrame(1,2);
t718 = sin(t743);
t849 = t664 * t718;
t722 = cos(t743);
t848 = t664 * t722;
t798 = t750 * mrSges(3,1) - mrSges(3,2) * t744;
t668 = t798 * t734 - (t744 * mrSges(3,1) + t750 * mrSges(3,2)) * t830;
t847 = t668 * t772;
t797 = t752 * mrSges(3,1) - mrSges(3,2) * t746;
t669 = t797 * t734 - (t746 * mrSges(3,1) + t752 * mrSges(3,2)) * t828;
t846 = t669 * t772;
t796 = t754 * mrSges(3,1) - mrSges(3,2) * t748;
t670 = t796 * t734 - (t748 * mrSges(3,1) + t754 * mrSges(3,2)) * t826;
t845 = t670 * t772;
t708 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t844 = ((t799 + t869) * t738 + t736 * t708) * t732;
t843 = ((t798 + t869) * t751 + t745 * t708) * t732;
t842 = ((t797 + t869) * t753 + t747 * t708) * t732;
t841 = ((t796 + t869) * t755 + t749 * t708) * t732;
t710 = -mrSges(3,2) * pkin(6) + Ifges(3,6);
t711 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t681 = t710 * t737 - t735 * t711;
t840 = t681 * t772;
t688 = t710 * t750 - t744 * t711;
t839 = t688 * t772;
t689 = t710 * t752 - t746 * t711;
t838 = t689 * t772;
t690 = t710 * t754 - t748 * t711;
t837 = t690 * t772;
t834 = t732 * t733;
t833 = t732 * t735;
t831 = t732 * t744;
t829 = t732 * t746;
t827 = t732 * t748;
t824 = t733 * t737;
t823 = t733 * t750;
t822 = t733 * t752;
t821 = t733 * t754;
t808 = -0.2e1 * pkin(2) * mrSges(3,2);
t807 = t715 * t867;
t806 = t719 * t867;
t805 = t716 * t862;
t804 = t720 * t862;
t803 = t717 * t860;
t802 = t721 * t860;
t801 = t718 * t858;
t800 = t722 * t858;
t761 = xP(4);
t723 = sin(t761);
t724 = cos(t761);
t762 = mrSges(4,2);
t763 = mrSges(4,1);
t795 = -t723 * t762 + t724 * t763;
t794 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t793 = -t723 * t763 - t724 * t762;
t792 = pkin(3) * t833 - t699 * t734;
t791 = pkin(3) * t831 - t701 * t734;
t790 = pkin(3) * t829 - t702 * t734;
t789 = pkin(3) * t827 - t703 * t734;
t788 = Ifges(3,3) * t866 + t660 * t681;
t787 = Ifges(3,3) * t861 + t662 * t688;
t786 = Ifges(3,3) * t859 + t663 * t689;
t785 = Ifges(3,3) * t857 + t664 * t690;
t739 = -Ifges(3,1) + Ifges(3,2);
t756 = mrSges(3,1) * pkin(2);
t655 = t739 * t725 + 0.2e1 * (Ifges(3,4) * t735 + t756) * t737 + t735 * t808 + t794;
t784 = t640 * t840 + t655 * t660;
t656 = t739 * t728 + 0.2e1 * (Ifges(3,4) * t744 + t756) * t750 + t744 * t808 + t794;
t783 = t644 * t839 + t656 * t662;
t657 = t739 * t729 + 0.2e1 * (Ifges(3,4) * t746 + t756) * t752 + t746 * t808 + t794;
t782 = t645 * t838 + t657 * t663;
t658 = t739 * t730 + 0.2e1 * (Ifges(3,4) * t748 + t756) * t754 + t748 * t808 + t794;
t781 = t646 * t837 + t658 * t664;
t764 = koppelP(4,2);
t768 = koppelP(4,1);
t691 = -t723 * t768 - t724 * t764;
t695 = -t723 * t764 + t724 * t768;
t780 = t647 * (-t691 * t719 + t695 * t715);
t765 = koppelP(3,2);
t769 = koppelP(3,1);
t692 = -t723 * t769 - t724 * t765;
t696 = -t723 * t765 + t724 * t769;
t779 = t648 * (-t692 * t720 + t696 * t716);
t766 = koppelP(2,2);
t770 = koppelP(2,1);
t693 = -t723 * t770 - t724 * t766;
t697 = -t723 * t766 + t724 * t770;
t778 = t649 * (-t693 * t721 + t697 * t717);
t767 = koppelP(1,2);
t771 = koppelP(1,1);
t694 = -t723 * t771 - t724 * t767;
t698 = -t723 * t767 + t724 * t771;
t777 = t650 * (-t694 * t722 + t698 * t718);
t776 = t640 * t856 + t660 * t844;
t775 = t644 * t847 + t662 * t843;
t774 = t645 * t846 + t663 * t842;
t773 = t646 * t845 + t664 * t841;
t726 = m(1) + m(2) + m(3);
t687 = t731 * t755 + t733 * t816;
t686 = t731 * t753 + t733 * t817;
t685 = t731 * t751 + t733 * t818;
t680 = t731 * t738 + t733 * t820;
t667 = -t748 * t687 - t732 * t821;
t666 = -t746 * t686 - t732 * t822;
t665 = -t744 * t685 - t732 * t823;
t661 = -t735 * t680 - t732 * t824;
t654 = -t706 * t731 + t789 * t733;
t653 = -t705 * t731 + t790 * t733;
t652 = -t704 * t731 + t791 * t733;
t651 = -t700 * t731 + t792 * t733;
t638 = -t684 * t875 + t706 * t821 + (pkin(2) * t827 + t789 * t754) * t731;
t637 = -t683 * t876 + t705 * t822 + (pkin(2) * t829 + t790 * t752) * t731;
t636 = -t682 * t877 + t704 * t823 + (pkin(2) * t831 + t791 * t750) * t731;
t635 = -t679 * t878 + t700 * t824 + (pkin(2) * t833 + t792 * t737) * t731;
t634 = -(t687 * t718 - t722 * t826) * t875 + (t654 * t718 + t722 * t678) * t754 + (t718 * t834 + t734 * t722) * t879;
t633 = -(t686 * t717 - t721 * t828) * t876 + (t653 * t717 + t721 * t677) * t752 + (t717 * t834 + t734 * t721) * t880;
t632 = -(t685 * t716 - t720 * t830) * t877 + (t652 * t716 + t720 * t676) * t750 + (t716 * t834 + t734 * t720) * t881;
t631 = (t687 * t722 + t718 * t826) * t875 + (-t654 * t722 + t718 * t678) * t754 + (t734 * t718 - t722 * t834) * t879;
t630 = (t686 * t721 + t717 * t828) * t876 + (-t653 * t721 + t717 * t677) * t752 + (t734 * t717 - t721 * t834) * t880;
t629 = (t685 * t720 + t716 * t830) * t877 + (-t652 * t720 + t716 * t676) * t750 + (t734 * t716 - t720 * t834) * t881;
t628 = -(t680 * t715 - t719 * t832) * t878 + (t651 * t715 + t719 * t675) * t737 + (t715 * t834 + t734 * t719) * t882;
t627 = (t680 * t719 + t715 * t832) * t878 + (-t651 * t719 + t715 * t675) * t737 + (t734 * t715 - t719 * t834) * t882;
t626 = t664 * t777;
t625 = t663 * t778;
t624 = t662 * t779;
t623 = t660 * t780;
t622 = t777 * t857;
t621 = t778 * t859;
t620 = t779 * t861;
t619 = t780 * t866;
t618 = (t638 * t670 + t643 * t870 + t667 * t690) * t650;
t617 = (t637 * t669 + t642 * t870 + t666 * t689) * t649;
t616 = (t636 * t668 + t641 * t870 + t665 * t688) * t648;
t615 = (t635 * t659 + t639 * t870 + t661 * t681) * t647;
t614 = (t638 * t726 + t643 * t845 + t667 * t841) * t650;
t613 = (t637 * t726 + t642 * t846 + t666 * t842) * t649;
t612 = (t636 * t726 + t641 * t847 + t665 * t843) * t648;
t611 = (t635 * t726 + t639 * t856 + t661 * t844) * t647;
t610 = (t638 * t841 + t643 * t837 + t658 * t667) * t650;
t609 = (t637 * t842 + t642 * t838 + t657 * t666) * t649;
t608 = (t636 * t843 + t641 * t839 + t656 * t665) * t648;
t607 = (t635 * t844 + t639 * t840 + t655 * t661) * t647;
t606 = (t631 * t694 + t634 * t698) * t650;
t605 = (t630 * t693 + t633 * t697) * t649;
t604 = (t629 * t692 + t632 * t696) * t648;
t603 = (t627 * t691 + t628 * t695) * t647;
t602 = (t634 * t670 + t785 * t718) * t650;
t601 = (t633 * t669 + t786 * t717) * t649;
t600 = (t632 * t668 + t787 * t716) * t648;
t599 = (t631 * t670 - t785 * t722) * t650;
t598 = (t630 * t669 - t786 * t721) * t649;
t597 = (t629 * t668 - t787 * t720) * t648;
t596 = (t628 * t659 + t788 * t715) * t647;
t595 = (t627 * t659 - t788 * t719) * t647;
t594 = (t634 * t726 + t773 * t718) * t650;
t593 = (t633 * t726 + t774 * t717) * t649;
t592 = (t632 * t726 + t775 * t716) * t648;
t591 = (t631 * t726 - t773 * t722) * t650;
t590 = (t630 * t726 - t774 * t721) * t649;
t589 = (t629 * t726 - t775 * t720) * t648;
t588 = (t628 * t726 + t776 * t715) * t647;
t587 = (t627 * t726 - t776 * t719) * t647;
t586 = (t634 * t841 + t781 * t718) * t650;
t585 = (t633 * t842 + t782 * t717) * t649;
t584 = (t632 * t843 + t783 * t716) * t648;
t583 = (t631 * t841 - t781 * t722) * t650;
t582 = (t630 * t842 - t782 * t721) * t649;
t581 = (t629 * t843 - t783 * t720) * t648;
t580 = (t628 * t844 + t784 * t715) * t647;
t579 = (t627 * t844 - t784 * t719) * t647;
t578 = t622 * Ifges(3,3) + t606 * t670 + t626 * t690;
t577 = t621 * Ifges(3,3) + t605 * t669 + t625 * t689;
t576 = t620 * Ifges(3,3) + t604 * t668 + t624 * t688;
t575 = t606 * t726 + t622 * t670 + t626 * t841;
t574 = t605 * t726 + t621 * t669 + t625 * t842;
t573 = t604 * t726 + t620 * t668 + t624 * t843;
t572 = t619 * Ifges(3,3) + t603 * t659 + t623 * t681;
t571 = t603 * t726 + t619 * t659 + t623 * t844;
t570 = t606 * t841 + t622 * t690 + t626 * t658;
t569 = t605 * t842 + t621 * t689 + t625 * t657;
t568 = t604 * t843 + t620 * t688 + t624 * t656;
t567 = t603 * t844 + t619 * t681 + t623 * t655;
t1 = [m(4) + (-t583 * t848 + t591 * t631) * t650 + (-t582 * t850 + t590 * t630) * t649 + (-t581 * t852 + t589 * t629) * t648 + (-t579 * t854 + t587 * t627) * t647 + (-t595 * t806 - t597 * t804 - t598 * t802 - t599 * t800) * t772, (t583 * t849 + t591 * t634) * t650 + (t582 * t851 + t590 * t633) * t649 + (t581 * t853 + t589 * t632) * t648 + (t579 * t855 + t587 * t628) * t647 + (t595 * t807 + t597 * t805 + t598 * t803 + t599 * t801) * t772, (t583 * t667 + t591 * t638) * t650 + (t582 * t666 + t590 * t637) * t649 + (t581 * t665 + t589 * t636) * t648 + (t579 * t661 + t587 * t635) * t647 + (t595 * t868 + t597 * t865 + t598 * t864 + t599 * t863) * t772, t579 * t623 + t581 * t624 + t582 * t625 + t583 * t626 + t587 * t603 + t589 * t604 + t590 * t605 + t591 * t606 + t595 * t619 + t597 * t620 + t598 * t621 + t599 * t622 + t793; (-t586 * t848 + t594 * t631) * t650 + (-t585 * t850 + t593 * t630) * t649 + (-t584 * t852 + t592 * t629) * t648 + (-t580 * t854 + t588 * t627) * t647 + (-t596 * t806 - t600 * t804 - t601 * t802 - t602 * t800) * t772, m(4) + (t586 * t849 + t594 * t634) * t650 + (t585 * t851 + t593 * t633) * t649 + (t584 * t853 + t592 * t632) * t648 + (t580 * t855 + t588 * t628) * t647 + (t596 * t807 + t600 * t805 + t601 * t803 + t602 * t801) * t772, (t586 * t667 + t594 * t638) * t650 + (t585 * t666 + t593 * t637) * t649 + (t584 * t665 + t592 * t636) * t648 + (t580 * t661 + t588 * t635) * t647 + (t596 * t868 + t600 * t865 + t601 * t864 + t602 * t863) * t772, t580 * t623 + t584 * t624 + t585 * t625 + t586 * t626 + t588 * t603 + t592 * t604 + t593 * t605 + t594 * t606 + t596 * t619 + t600 * t620 + t601 * t621 + t602 * t622 + t795; (-t610 * t848 + t614 * t631) * t650 + (-t609 * t850 + t613 * t630) * t649 + (-t608 * t852 + t612 * t629) * t648 + (-t607 * t854 + t611 * t627) * t647 + (-t615 * t806 - t616 * t804 - t617 * t802 - t618 * t800) * t772, (t610 * t849 + t614 * t634) * t650 + (t609 * t851 + t613 * t633) * t649 + (t608 * t853 + t612 * t632) * t648 + (t607 * t855 + t611 * t628) * t647 + (t615 * t807 + t616 * t805 + t617 * t803 + t618 * t801) * t772, m(4) + (t610 * t667 + t614 * t638) * t650 + (t609 * t666 + t613 * t637) * t649 + (t608 * t665 + t612 * t636) * t648 + (t607 * t661 + t611 * t635) * t647 + (t615 * t868 + t616 * t865 + t617 * t864 + t618 * t863) * t772, t611 * t603 + t612 * t604 + t613 * t605 + t614 * t606 + t607 * t623 + t608 * t624 + t609 * t625 + t610 * t626 + t615 * t619 + t616 * t620 + t617 * t621 + t618 * t622; (-t570 * t848 + t575 * t631) * t650 + (-t569 * t850 + t574 * t630) * t649 + (-t568 * t852 + t573 * t629) * t648 + (-t567 * t854 + t571 * t627) * t647 + (-t572 * t806 - t576 * t804 - t577 * t802 - t578 * t800) * t772 + t793, (t570 * t849 + t575 * t634) * t650 + (t569 * t851 + t574 * t633) * t649 + (t568 * t853 + t573 * t632) * t648 + (t567 * t855 + t571 * t628) * t647 + (t572 * t807 + t576 * t805 + t577 * t803 + t578 * t801) * t772 + t795, (t570 * t667 + t575 * t638) * t650 + (t569 * t666 + t574 * t637) * t649 + (t568 * t665 + t573 * t636) * t648 + (t567 * t661 + t571 * t635) * t647 + (t572 * t868 + t576 * t865 + t577 * t864 + t578 * t863) * t772, t567 * t623 + t568 * t624 + t569 * t625 + t570 * t626 + t571 * t603 + t572 * t619 + t573 * t604 + t574 * t605 + t575 * t606 + t576 * t620 + t577 * t621 + t578 * t622 + Ifges(4,3);];
MX  = t1;
