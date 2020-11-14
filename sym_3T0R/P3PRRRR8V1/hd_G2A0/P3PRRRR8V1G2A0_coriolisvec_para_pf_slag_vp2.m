% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 17:03
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:02:52
% EndTime: 2020-08-06 17:02:55
% DurationCPUTime: 3.44s
% Computational Cost: add. (15144->217), mult. (48621->461), div. (6885->7), fcn. (54228->22), ass. (0->192)
t702 = cos(qJ(3,1));
t703 = cos(qJ(2,1));
t697 = sin(qJ(2,1));
t749 = t697 * t702;
t654 = pkin(2) * t749 - t703 * pkin(5);
t684 = sin(pkin(3));
t696 = sin(qJ(3,1));
t686 = cos(pkin(3));
t793 = pkin(2) * t686;
t645 = t654 * t684 + t696 * t793;
t806 = 0.1e1 / t645;
t774 = t806 / t702;
t700 = cos(qJ(3,2));
t701 = cos(qJ(2,2));
t695 = sin(qJ(2,2));
t752 = t695 * t700;
t653 = pkin(2) * t752 - t701 * pkin(5);
t694 = sin(qJ(3,2));
t644 = t653 * t684 + t694 * t793;
t807 = 0.1e1 / t644;
t775 = t807 / t700;
t698 = cos(qJ(3,3));
t699 = cos(qJ(2,3));
t693 = sin(qJ(2,3));
t755 = t693 * t698;
t652 = pkin(2) * t755 - t699 * pkin(5);
t692 = sin(qJ(3,3));
t643 = t652 * t684 + t692 * t793;
t808 = 0.1e1 / t643;
t776 = t808 / t698;
t736 = t698 * mrSges(3,1) - t692 * mrSges(3,2);
t735 = t700 * mrSges(3,1) - t694 * mrSges(3,2);
t734 = t702 * mrSges(3,1) - t696 * mrSges(3,2);
t809 = -2 * Ifges(3,4);
t677 = t698 ^ 2;
t805 = 0.2e1 * t677;
t679 = t700 ^ 2;
t804 = 0.2e1 * t679;
t681 = t702 ^ 2;
t803 = 0.2e1 * t681;
t802 = Ifges(3,5) / 0.2e1;
t801 = -Ifges(3,6) / 0.2e1;
t688 = mrSges(2,2) - mrSges(3,3);
t800 = t688 / 0.2e1;
t683 = sin(pkin(6));
t685 = cos(pkin(6));
t761 = t686 * t699;
t764 = t686 * t693;
t792 = pkin(2) * t698;
t622 = -(-t683 * t693 + t685 * t761) * t792 - pkin(5) * (t683 * t699 + t685 * t764);
t625 = (t683 * t761 + t685 * t693) * t792 + (t683 * t764 - t685 * t699) * pkin(5);
t704 = xDP(3);
t709 = 0.1e1 / pkin(2);
t689 = legFrame(3,2);
t670 = sin(t689);
t673 = cos(t689);
t705 = xDP(2);
t706 = xDP(1);
t724 = t670 * t705 - t673 * t706;
t601 = (t622 * t704 + t724 * t625) * t709 * t776;
t799 = pkin(2) * t601;
t760 = t686 * t701;
t763 = t686 * t695;
t791 = pkin(2) * t700;
t623 = -(-t683 * t695 + t685 * t760) * t791 - pkin(5) * (t683 * t701 + t685 * t763);
t626 = (t683 * t760 + t685 * t695) * t791 + (t683 * t763 - t685 * t701) * pkin(5);
t690 = legFrame(2,2);
t671 = sin(t690);
t674 = cos(t690);
t723 = t671 * t705 - t674 * t706;
t602 = (t623 * t704 + t723 * t626) * t709 * t775;
t798 = pkin(2) * t602;
t759 = t686 * t703;
t762 = t686 * t697;
t790 = pkin(2) * t702;
t624 = -(-t683 * t697 + t685 * t759) * t790 - pkin(5) * (t683 * t703 + t685 * t762);
t627 = (t683 * t759 + t685 * t697) * t790 + (t683 * t762 - t685 * t703) * pkin(5);
t691 = legFrame(1,2);
t672 = sin(t691);
t675 = cos(t691);
t722 = t672 * t705 - t675 * t706;
t603 = (t624 * t704 + t722 * t627) * t709 * t774;
t797 = pkin(2) * t603;
t796 = pkin(2) * t677;
t795 = pkin(2) * t679;
t794 = pkin(2) * t681;
t789 = -Ifges(3,1) - Ifges(2,3);
t712 = t684 * t698 + t692 * t764;
t756 = t692 * t699;
t628 = t712 * t683 - t685 * t756;
t631 = -t683 * t756 - t712 * t685;
t607 = (t724 * t628 + t631 * t704) * t776;
t707 = pkin(5) ^ 2;
t708 = pkin(2) ^ 2;
t782 = t601 * t692;
t748 = pkin(2) * t782;
t785 = (-pkin(5) * t748 + (t677 * t708 + t707) * t607) * t607;
t711 = t684 * t700 + t694 * t763;
t753 = t694 * t701;
t629 = t711 * t683 - t685 * t753;
t632 = -t683 * t753 - t711 * t685;
t608 = (t723 * t629 + t632 * t704) * t775;
t781 = t602 * t694;
t747 = pkin(2) * t781;
t784 = (-pkin(5) * t747 + (t679 * t708 + t707) * t608) * t608;
t710 = t684 * t702 + t696 * t762;
t750 = t696 * t703;
t630 = t710 * t683 - t685 * t750;
t633 = -t683 * t750 - t710 * t685;
t609 = (t722 * t630 + t633 * t704) * t774;
t780 = t603 * t696;
t746 = pkin(2) * t780;
t783 = (-pkin(5) * t746 + (t681 * t708 + t707) * t609) * t609;
t779 = t607 * t692;
t778 = t608 * t694;
t777 = t609 * t696;
t773 = ((mrSges(2,1) + t736) * t699 - t693 * t688) * t684;
t772 = ((mrSges(2,1) + t735) * t701 - t695 * t688) * t684;
t771 = ((mrSges(2,1) + t734) * t703 - t697 * t688) * t684;
t770 = t684 * t692;
t769 = t684 * t694;
t768 = t684 * t696;
t767 = t684 * t699;
t766 = t684 * t701;
t765 = t684 * t703;
t758 = t686 * t709;
t757 = t692 * t698;
t754 = t694 * t700;
t751 = t696 * t702;
t745 = pkin(5) * t779;
t744 = pkin(5) * t778;
t743 = pkin(5) * t777;
t742 = t673 * t776;
t741 = t674 * t775;
t740 = t675 * t774;
t739 = t684 * t755;
t738 = t684 * t752;
t737 = t684 * t749;
t592 = t745 - t799;
t574 = (((t686 * t601 + t607 * t767) * t796 - (-pkin(5) * t607 + t748) * t739 + t686 * t592) * t607 - (-t601 * t767 + (-t677 * t686 + t692 * t739 + t686) * t607) * t799) * t776;
t580 = (t758 * t785 + (-t601 * t652 * t770 + t686 * (t601 * t796 - t745)) * t601) * t776;
t583 = (t592 * t799 - t785) * t808;
t598 = t601 ^ 2;
t604 = t607 ^ 2;
t658 = t692 * mrSges(3,1) + t698 * mrSges(3,2);
t634 = t693 * t658 * t684 - t736 * t686;
t676 = -m(1) - m(2) - m(3);
t733 = t808 * (-t574 * t773 + t634 * t580 + t676 * t583 + ((-t604 * mrSges(2,1) - t736 * (t604 + t598)) * t693 - 0.2e1 * (t658 * t601 + t607 * t800) * t607 * t699) * t684 - t686 * t598 * t658);
t594 = -t743 + t797;
t576 = (((t686 * t603 + t609 * t765) * t794 - (-pkin(5) * t609 + t746) * t737 - t686 * t594) * t609 - (-t603 * t765 + (-t681 * t686 + t696 * t737 + t686) * t609) * t797) * t774;
t582 = (t758 * t783 + (-t603 * t654 * t768 + t686 * (t603 * t794 - t743)) * t603) * t774;
t585 = (-t594 * t797 - t783) * t806;
t600 = t603 ^ 2;
t606 = t609 ^ 2;
t660 = t696 * mrSges(3,1) + t702 * mrSges(3,2);
t636 = t697 * t660 * t684 - t734 * t686;
t732 = t806 * (-t576 * t771 + t636 * t582 + t676 * t585 + ((-t606 * mrSges(2,1) - t734 * (t606 + t600)) * t697 - 0.2e1 * (t660 * t603 + t609 * t800) * t609 * t703) * t684 - t686 * t600 * t660);
t593 = t744 - t798;
t575 = (((t686 * t602 + t608 * t766) * t795 - (-pkin(5) * t608 + t747) * t738 + t686 * t593) * t608 + (t602 * t766 + (t679 * t686 - t694 * t738 - t686) * t608) * t798) * t775;
t581 = (t758 * t784 + (-t602 * t653 * t769 + t686 * (t602 * t795 - t744)) * t602) * t775;
t584 = (t593 * t798 - t784) * t807;
t599 = t602 ^ 2;
t605 = t608 ^ 2;
t659 = t694 * mrSges(3,1) + t700 * mrSges(3,2);
t635 = t695 * t659 * t684 - t735 * t686;
t731 = (-t575 * t772 + t635 * t581 + t676 * t584 + ((-t605 * mrSges(2,1) - t735 * (t605 + t599)) * t695 - 0.2e1 * (t659 * t602 + t608 * t800) * t608 * t701) * t684 - t686 * t599 * t659) * t807;
t661 = -Ifges(3,5) * t692 - Ifges(3,6) * t698;
t687 = Ifges(3,1) - Ifges(3,2);
t730 = 0.2e1 * ((t601 * t802 + t687 * t779) * t698 + t782 * t801 + (t805 - 0.1e1) * t607 * Ifges(3,4)) * t601 - t583 * t773 + (t687 * t677 + t757 * t809 + t789) * t574 + t661 * t580;
t662 = -Ifges(3,5) * t694 - Ifges(3,6) * t700;
t729 = 0.2e1 * ((t602 * t802 + t687 * t778) * t700 + t781 * t801 + (t804 - 0.1e1) * t608 * Ifges(3,4)) * t602 - t584 * t772 + (t687 * t679 + t754 * t809 + t789) * t575 + t662 * t581;
t663 = -Ifges(3,5) * t696 - Ifges(3,6) * t702;
t728 = 0.2e1 * ((t603 * t802 + t687 * t777) * t702 + t780 * t801 + (t803 - 0.1e1) * t609 * Ifges(3,4)) * t603 - t585 * t771 + (t687 * t681 + t751 * t809 + t789) * t576 + t663 * t582;
t727 = Ifges(3,3) * t580 - t661 * t574 - t634 * t583 + t604 * (Ifges(3,4) * t805 + t687 * t757 - Ifges(3,4));
t726 = Ifges(3,3) * t581 - t662 * t575 - t635 * t584 + t605 * (Ifges(3,4) * t804 + t687 * t754 - Ifges(3,4));
t725 = Ifges(3,3) * t582 - t663 * t576 - t636 * t585 + t606 * (Ifges(3,4) * t803 + t687 * t751 - Ifges(3,4));
t721 = pkin(2) * t770 - t652 * t686;
t720 = pkin(2) * t769 - t653 * t686;
t719 = pkin(2) * t768 - t654 * t686;
t718 = t730 * t776;
t717 = t727 * t776;
t716 = t729 * t775;
t715 = t726 * t775;
t714 = t728 * t774;
t713 = t725 * t774;
t657 = pkin(5) * t697 + t703 * t790;
t656 = pkin(5) * t695 + t701 * t791;
t655 = pkin(5) * t693 + t699 * t792;
t618 = -t683 * t657 + t719 * t685;
t617 = -t683 * t656 + t720 * t685;
t616 = -t683 * t655 + t721 * t685;
t1 = [-t728 * t630 * t740 - t729 * t629 * t741 - t730 * t628 * t742 + (-t618 * t675 + t645 * t672) * t732 + (-t617 * t674 + t644 * t671) * t731 + (-t616 * t673 + t643 * t670) * t733 + (t727 * t625 * t742 + t726 * t626 * t741 + t725 * t627 * t740) * t709; t672 * t630 * t714 + t671 * t629 * t716 + t670 * t628 * t718 + (t618 * t672 + t645 * t675) * t732 + (t617 * t671 + t644 * t674) * t731 + (t616 * t670 + t643 * t673) * t733 + (-t670 * t625 * t717 - t671 * t626 * t715 - t672 * t627 * t713) * t709; t633 * t714 + t632 * t716 + t631 * t718 + (t685 * t657 + t719 * t683) * t732 + (t685 * t656 + t720 * t683) * t731 + (t685 * t655 + t721 * t683) * t733 + (-t622 * t717 - t623 * t715 - t624 * t713) * t709;];
taucX  = t1;
