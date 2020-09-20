% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G1A0
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:49:57
% EndTime: 2020-08-06 16:50:01
% DurationCPUTime: 3.23s
% Computational Cost: add. (13098->231), mult. (34998->476), div. (3729->8), fcn. (42558->22), ass. (0->212)
t696 = cos(qJ(3,1));
t676 = 0.1e1 / t696;
t697 = cos(qJ(2,1));
t691 = sin(qJ(2,1));
t738 = t691 * t696;
t649 = pkin(2) * t738 - t697 * pkin(5);
t678 = sin(pkin(3));
t690 = sin(qJ(3,1));
t680 = cos(pkin(3));
t789 = pkin(2) * t680;
t712 = 0.1e1 / (t649 * t678 + t690 * t789);
t767 = t712 * t676;
t694 = cos(qJ(3,2));
t674 = 0.1e1 / t694;
t695 = cos(qJ(2,2));
t689 = sin(qJ(2,2));
t742 = t689 * t694;
t648 = pkin(2) * t742 - t695 * pkin(5);
t688 = sin(qJ(3,2));
t713 = 0.1e1 / (t648 * t678 + t688 * t789);
t768 = t713 * t674;
t681 = legFrame(3,3);
t661 = sin(t681);
t664 = cos(t681);
t677 = sin(pkin(6));
t679 = cos(pkin(6));
t632 = -t677 * t661 + t664 * t679;
t693 = cos(qJ(2,3));
t687 = sin(qJ(2,3));
t753 = t680 * t687;
t638 = t677 * t693 + t679 * t753;
t641 = -t677 * t753 + t679 * t693;
t686 = sin(qJ(3,3));
t692 = cos(qJ(3,3));
t760 = t678 * t692;
t606 = (-t638 * t664 - t661 * t641) * t686 - t632 * t760;
t635 = t679 * t661 + t664 * t677;
t609 = (-t661 * t638 + t641 * t664) * t686 - t635 * t760;
t746 = t687 * t692;
t647 = pkin(2) * t746 - t693 * pkin(5);
t754 = t680 * t686;
t626 = 0.1e1 / (pkin(2) * t754 + t647 * t678);
t672 = 0.1e1 / t692;
t698 = xDP(2);
t699 = xDP(1);
t591 = (t606 * t699 + t609 * t698) * t672 * t626;
t671 = t692 ^ 2;
t700 = pkin(5) ^ 2;
t701 = pkin(2) ^ 2;
t736 = t693 * t635;
t737 = t693 * t632;
t747 = t687 * t635;
t748 = t687 * t632;
t788 = pkin(2) * t692;
t600 = -(t680 * t737 - t747) * t788 - pkin(5) * (t680 * t748 + t736);
t601 = -(t680 * t736 + t748) * t788 - (t680 * t747 - t737) * pkin(5);
t702 = 0.1e1 / pkin(2);
t719 = t678 * t746;
t759 = t678 * t693;
t802 = 0.1e1 / (-pkin(5) * t759 + (t719 + t754) * pkin(2));
t769 = t802 * t672;
t585 = (t600 * t699 + t601 * t698) * t702 * t769;
t776 = t585 * t686;
t725 = pkin(2) * t776;
t573 = -pkin(5) * t725 + (t671 * t701 + t700) * t591;
t772 = t591 * t686;
t722 = pkin(5) * t772;
t750 = t680 * t702;
t763 = t678 * t686;
t773 = t591 * t802;
t792 = pkin(2) * t671;
t561 = (t573 * t750 * t773 + (-t585 * t647 * t763 + t680 * (t585 * t792 - t722)) * t626 * t585) * t672;
t795 = pkin(2) * t585;
t576 = t722 - t795;
t567 = (t573 * t591 - t576 * t795) * t802;
t582 = t585 ^ 2;
t588 = t591 ^ 2;
t653 = t686 * mrSges(3,1) + t692 * mrSges(3,2);
t716 = t692 * mrSges(3,1) - t686 * mrSges(3,2);
t618 = t653 * t678 * t687 - t680 * t716;
t670 = -m(1) - m(2) - m(3);
t685 = mrSges(2,2) - mrSges(3,3);
t796 = t685 / 0.2e1;
t806 = ((-t588 * mrSges(2,1) - t716 * (t588 + t582)) * t687 - 0.2e1 * t591 * (t585 * t653 + t591 * t796) * t693) * t678 + t618 * t561 - t670 * t567;
t682 = legFrame(2,3);
t662 = sin(t682);
t665 = cos(t682);
t636 = t679 * t662 + t665 * t677;
t734 = t695 * t636;
t633 = -t677 * t662 + t665 * t679;
t735 = t695 * t633;
t743 = t689 * t636;
t744 = t689 * t633;
t787 = pkin(2) * t694;
t602 = -(t680 * t735 - t743) * t787 - pkin(5) * (t680 * t744 + t734);
t603 = -(t680 * t734 + t744) * t787 - (t680 * t743 - t735) * pkin(5);
t586 = (t602 * t699 + t603 * t698) * t702 * t768;
t752 = t680 * t689;
t639 = t677 * t695 + t679 * t752;
t642 = -t677 * t752 + t679 * t695;
t758 = t678 * t694;
t607 = (-t639 * t665 - t662 * t642) * t688 - t633 * t758;
t610 = (-t662 * t639 + t642 * t665) * t688 - t636 * t758;
t592 = (t607 * t699 + t610 * t698) * t768;
t771 = t592 * t688;
t721 = pkin(5) * t771;
t762 = t678 * t688;
t673 = t694 ^ 2;
t775 = t586 * t688;
t724 = pkin(2) * t775;
t781 = (-pkin(5) * t724 + (t673 * t701 + t700) * t592) * t592;
t791 = pkin(2) * t673;
t562 = (t750 * t781 + (-t586 * t648 * t762 + t680 * (t586 * t791 - t721)) * t586) * t768;
t794 = pkin(2) * t586;
t577 = t721 - t794;
t568 = (t577 * t794 - t781) * t713;
t583 = t586 ^ 2;
t589 = t592 ^ 2;
t654 = t688 * mrSges(3,1) + t694 * mrSges(3,2);
t715 = t694 * mrSges(3,1) - t688 * mrSges(3,2);
t619 = t654 * t678 * t689 - t680 * t715;
t805 = ((-t589 * mrSges(2,1) - t715 * (t589 + t583)) * t689 - 0.2e1 * t592 * (t586 * t654 + t592 * t796) * t695) * t678 + t619 * t562 + t670 * t568;
t683 = legFrame(1,3);
t663 = sin(t683);
t666 = cos(t683);
t637 = t679 * t663 + t666 * t677;
t732 = t697 * t637;
t634 = -t677 * t663 + t666 * t679;
t733 = t697 * t634;
t739 = t691 * t637;
t740 = t691 * t634;
t786 = pkin(2) * t696;
t604 = -(t680 * t733 - t739) * t786 - pkin(5) * (t680 * t740 + t732);
t605 = -(t680 * t732 + t740) * t786 - (t680 * t739 - t733) * pkin(5);
t587 = (t604 * t699 + t605 * t698) * t702 * t767;
t751 = t680 * t691;
t640 = t677 * t697 + t679 * t751;
t643 = -t677 * t751 + t679 * t697;
t756 = t678 * t696;
t608 = (-t640 * t666 - t663 * t643) * t690 - t634 * t756;
t611 = (-t663 * t640 + t643 * t666) * t690 - t637 * t756;
t593 = (t608 * t699 + t611 * t698) * t767;
t770 = t593 * t690;
t720 = pkin(5) * t770;
t761 = t678 * t690;
t675 = t696 ^ 2;
t774 = t587 * t690;
t723 = pkin(2) * t774;
t780 = (-pkin(5) * t723 + (t675 * t701 + t700) * t593) * t593;
t790 = pkin(2) * t675;
t563 = (t750 * t780 + (-t587 * t649 * t761 + t680 * (t587 * t790 - t720)) * t587) * t767;
t793 = pkin(2) * t587;
t578 = t720 - t793;
t569 = (t578 * t793 - t780) * t712;
t584 = t587 ^ 2;
t590 = t593 ^ 2;
t655 = t690 * mrSges(3,1) + t696 * mrSges(3,2);
t714 = t696 * mrSges(3,1) - t690 * mrSges(3,2);
t620 = t655 * t678 * t691 - t680 * t714;
t804 = ((-t590 * mrSges(2,1) - t714 * (t590 + t584)) * t691 - 0.2e1 * t593 * (t587 * t655 + t593 * t796) * t697) * t678 + t620 * t563 + t670 * t569;
t803 = -2 * Ifges(3,4);
t801 = 0.2e1 * t671;
t800 = 0.2e1 * t673;
t799 = 0.2e1 * t675;
t798 = Ifges(3,5) / 0.2e1;
t797 = -Ifges(3,6) / 0.2e1;
t785 = -Ifges(3,1) - Ifges(2,3);
t779 = t582 * t653;
t778 = t583 * t654;
t777 = t584 * t655;
t629 = (mrSges(2,1) + t716) * t693 - t687 * t685;
t766 = t629 * t678;
t630 = (mrSges(2,1) + t715) * t695 - t689 * t685;
t765 = t630 * t678;
t631 = (mrSges(2,1) + t714) * t697 - t691 * t685;
t764 = t631 * t678;
t757 = t678 * t695;
t755 = t678 * t697;
t749 = t686 * t692;
t745 = t688 * t694;
t741 = t690 * t696;
t549 = (((t680 * t585 + t591 * t759) * t792 - (-pkin(5) * t591 + t725) * t719 + t680 * t576) * t773 + (t585 * t759 + (t671 * t680 - t686 * t719 - t680) * t591) * t802 * t795) * t672;
t731 = -t549 * t766 - t680 * t779 + t806;
t718 = t678 * t742;
t550 = (((t680 * t586 + t592 * t757) * t791 - (-pkin(5) * t592 + t724) * t718 + t680 * t577) * t592 + (t586 * t757 + (t673 * t680 - t688 * t718 - t680) * t592) * t794) * t768;
t730 = -t550 * t765 - t680 * t778 + t805;
t717 = t678 * t738;
t551 = (((t680 * t587 + t593 * t755) * t790 - (-pkin(5) * t593 + t723) * t717 + t680 * t578) * t593 - (-t587 * t755 + (-t675 * t680 + t690 * t717 + t680) * t593) * t793) * t767;
t729 = -t551 * t764 - t680 * t777 + t804;
t656 = -Ifges(3,5) * t686 - Ifges(3,6) * t692;
t684 = Ifges(3,1) - Ifges(3,2);
t711 = (0.2e1 * ((t585 * t798 + t684 * t772) * t692 + t776 * t797 + (t801 - 0.1e1) * t591 * Ifges(3,4)) * t585 + t567 * t766 + (t684 * t671 + t749 * t803 + t785) * t549 + t656 * t561) * t672;
t657 = -Ifges(3,5) * t688 - Ifges(3,6) * t694;
t710 = (0.2e1 * ((t586 * t798 + t684 * t771) * t694 + t775 * t797 + (t800 - 0.1e1) * t592 * Ifges(3,4)) * t586 - t568 * t765 + (t684 * t673 + t745 * t803 + t785) * t550 + t657 * t562) * t674;
t658 = -Ifges(3,5) * t690 - Ifges(3,6) * t696;
t709 = (0.2e1 * ((t587 * t798 + t684 * t770) * t696 + t774 * t797 + (t799 - 0.1e1) * t593 * Ifges(3,4)) * t587 - t569 * t764 + (t684 * t675 + t741 * t803 + t785) * t551 + t658 * t563) * t676;
t708 = pkin(2) * t763 - t647 * t680;
t707 = pkin(2) * t762 - t648 * t680;
t706 = pkin(2) * t761 - t649 * t680;
t705 = (-Ifges(3,3) * t561 + t656 * t549 - t618 * t567 - t588 * (Ifges(3,4) * t801 + t684 * t749 - Ifges(3,4))) * t769;
t704 = (-Ifges(3,3) * t562 + t657 * t550 + t619 * t568 - t589 * (Ifges(3,4) * t800 + t684 * t745 - Ifges(3,4))) * t768;
t703 = (-Ifges(3,3) * t563 + t658 * t551 + t620 * t569 - t590 * (Ifges(3,4) * t799 + t684 * t741 - Ifges(3,4))) * t767;
t652 = pkin(5) * t691 + t697 * t786;
t651 = pkin(5) * t689 + t695 * t787;
t650 = pkin(5) * t687 + t693 * t788;
t617 = t677 * t652 - t679 * t706;
t616 = t677 * t651 - t679 * t707;
t615 = t677 * t650 - t679 * t708;
t614 = t652 * t679 + t677 * t706;
t613 = t651 * t679 + t677 * t707;
t612 = t650 * t679 + t677 * t708;
t1 = [(t608 * t709 + t729 * (t614 * t666 - t663 * t617)) * t712 + (t607 * t710 + t730 * (t613 * t665 - t662 * t616)) * t713 + (t606 * t711 + t731 * (t612 * t664 - t661 * t615)) * t626 + (t600 * t705 + t602 * t704 + t604 * t703) * t702; (t611 * t709 + t729 * (t614 * t663 + t617 * t666)) * t712 + (t610 * t710 + t730 * (t613 * t662 + t616 * t665)) * t713 + (t609 * t711 + t731 * (t612 * t661 + t615 * t664)) * t626 + (t601 * t705 + t603 * t704 + t605 * t703) * t702; (-t777 - t778 - t779) * t680 + (-t549 * t629 - t550 * t630 - t551 * t631) * t678 + t804 + t805 + t806;];
taucX  = t1;
