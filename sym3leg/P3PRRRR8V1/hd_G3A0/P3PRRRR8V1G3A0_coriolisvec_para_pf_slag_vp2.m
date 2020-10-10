% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRRR8V1G3A0
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:11
% EndTime: 2020-08-06 17:16:14
% DurationCPUTime: 3.44s
% Computational Cost: add. (15144->217), mult. (48621->461), div. (6885->7), fcn. (54228->22), ass. (0->192)
t699 = cos(qJ(3,1));
t700 = cos(qJ(2,1));
t694 = sin(qJ(2,1));
t746 = t694 * t699;
t651 = pkin(2) * t746 - pkin(5) * t700;
t681 = sin(pkin(3));
t693 = sin(qJ(3,1));
t683 = cos(pkin(3));
t790 = pkin(2) * t683;
t642 = t651 * t681 + t693 * t790;
t803 = 0.1e1 / t642;
t771 = t803 / t699;
t697 = cos(qJ(3,2));
t698 = cos(qJ(2,2));
t692 = sin(qJ(2,2));
t749 = t692 * t697;
t650 = pkin(2) * t749 - pkin(5) * t698;
t691 = sin(qJ(3,2));
t641 = t650 * t681 + t691 * t790;
t804 = 0.1e1 / t641;
t772 = t804 / t697;
t695 = cos(qJ(3,3));
t696 = cos(qJ(2,3));
t690 = sin(qJ(2,3));
t752 = t690 * t695;
t649 = pkin(2) * t752 - pkin(5) * t696;
t689 = sin(qJ(3,3));
t640 = t649 * t681 + t689 * t790;
t805 = 0.1e1 / t640;
t773 = t805 / t695;
t733 = t695 * mrSges(3,1) - mrSges(3,2) * t689;
t732 = t697 * mrSges(3,1) - mrSges(3,2) * t691;
t731 = t699 * mrSges(3,1) - mrSges(3,2) * t693;
t806 = -2 * Ifges(3,4);
t674 = t695 ^ 2;
t802 = 0.2e1 * t674;
t676 = t697 ^ 2;
t801 = 0.2e1 * t676;
t678 = t699 ^ 2;
t800 = 0.2e1 * t678;
t799 = Ifges(3,5) / 0.2e1;
t798 = -Ifges(3,6) / 0.2e1;
t685 = mrSges(2,2) - mrSges(3,3);
t797 = t685 / 0.2e1;
t680 = sin(pkin(6));
t682 = cos(pkin(6));
t758 = t683 * t696;
t761 = t683 * t690;
t789 = pkin(2) * t695;
t619 = (-t680 * t690 + t682 * t758) * t789 + pkin(5) * (t680 * t696 + t682 * t761);
t622 = (t680 * t758 + t682 * t690) * t789 + (t680 * t761 - t682 * t696) * pkin(5);
t701 = xDP(3);
t706 = 0.1e1 / pkin(2);
t686 = legFrame(3,2);
t667 = sin(t686);
t670 = cos(t686);
t702 = xDP(2);
t703 = xDP(1);
t721 = t667 * t702 - t670 * t703;
t598 = (t721 * t619 + t622 * t701) * t706 * t773;
t796 = pkin(2) * t598;
t757 = t683 * t698;
t760 = t683 * t692;
t788 = pkin(2) * t697;
t620 = (-t680 * t692 + t682 * t757) * t788 + pkin(5) * (t680 * t698 + t682 * t760);
t623 = (t680 * t757 + t682 * t692) * t788 + (t680 * t760 - t682 * t698) * pkin(5);
t687 = legFrame(2,2);
t668 = sin(t687);
t671 = cos(t687);
t720 = t668 * t702 - t671 * t703;
t599 = (t720 * t620 + t623 * t701) * t706 * t772;
t795 = pkin(2) * t599;
t756 = t683 * t700;
t759 = t683 * t694;
t787 = pkin(2) * t699;
t621 = (-t680 * t694 + t682 * t756) * t787 + pkin(5) * (t680 * t700 + t682 * t759);
t624 = (t680 * t756 + t682 * t694) * t787 + (t680 * t759 - t682 * t700) * pkin(5);
t688 = legFrame(1,2);
t669 = sin(t688);
t672 = cos(t688);
t719 = t669 * t702 - t672 * t703;
t600 = (t719 * t621 + t624 * t701) * t706 * t771;
t794 = pkin(2) * t600;
t793 = pkin(2) * t674;
t792 = pkin(2) * t676;
t791 = pkin(2) * t678;
t786 = -Ifges(3,1) - Ifges(2,3);
t709 = t681 * t695 + t689 * t761;
t753 = t689 * t696;
t625 = t709 * t680 - t682 * t753;
t628 = t680 * t753 + t709 * t682;
t604 = (t625 * t701 + t721 * t628) * t773;
t704 = pkin(5) ^ 2;
t705 = pkin(2) ^ 2;
t779 = t598 * t689;
t745 = pkin(2) * t779;
t782 = (-pkin(5) * t745 + (t674 * t705 + t704) * t604) * t604;
t708 = t681 * t697 + t691 * t760;
t750 = t691 * t698;
t626 = t708 * t680 - t682 * t750;
t629 = t680 * t750 + t708 * t682;
t605 = (t626 * t701 + t720 * t629) * t772;
t778 = t599 * t691;
t744 = pkin(2) * t778;
t781 = (-pkin(5) * t744 + (t676 * t705 + t704) * t605) * t605;
t707 = t681 * t699 + t693 * t759;
t747 = t693 * t700;
t627 = t707 * t680 - t682 * t747;
t630 = t680 * t747 + t707 * t682;
t606 = (t627 * t701 + t719 * t630) * t771;
t777 = t600 * t693;
t743 = pkin(2) * t777;
t780 = (-pkin(5) * t743 + (t678 * t705 + t704) * t606) * t606;
t776 = t604 * t689;
t775 = t605 * t691;
t774 = t606 * t693;
t770 = ((mrSges(2,1) + t733) * t696 - t690 * t685) * t681;
t769 = ((mrSges(2,1) + t732) * t698 - t692 * t685) * t681;
t768 = ((mrSges(2,1) + t731) * t700 - t694 * t685) * t681;
t767 = t681 * t689;
t766 = t681 * t691;
t765 = t681 * t693;
t764 = t681 * t696;
t763 = t681 * t698;
t762 = t681 * t700;
t755 = t683 * t706;
t754 = t689 * t695;
t751 = t691 * t697;
t748 = t693 * t699;
t742 = pkin(5) * t776;
t741 = pkin(5) * t775;
t740 = pkin(5) * t774;
t739 = t670 * t773;
t738 = t671 * t772;
t737 = t672 * t771;
t736 = t681 * t752;
t735 = t681 * t749;
t734 = t681 * t746;
t589 = t742 - t796;
t571 = (((t598 * t683 + t604 * t764) * t793 - (-pkin(5) * t604 + t745) * t736 + t683 * t589) * t604 - (-t598 * t764 + (-t674 * t683 + t689 * t736 + t683) * t604) * t796) * t773;
t577 = (t755 * t782 + (-t598 * t649 * t767 + t683 * (t598 * t793 - t742)) * t598) * t773;
t580 = (t589 * t796 - t782) * t805;
t595 = t598 ^ 2;
t601 = t604 ^ 2;
t655 = mrSges(3,1) * t689 + mrSges(3,2) * t695;
t631 = t681 * t655 * t690 - t733 * t683;
t673 = -m(1) - m(2) - m(3);
t730 = (-t571 * t770 + t577 * t631 + t580 * t673 + ((-mrSges(2,1) * t601 - t733 * (t601 + t595)) * t690 - 0.2e1 * t696 * t604 * (t655 * t598 + t604 * t797)) * t681 - t683 * t595 * t655) * t805;
t590 = t741 - t795;
t572 = (((t599 * t683 + t605 * t763) * t792 - (-pkin(5) * t605 + t744) * t735 + t683 * t590) * t605 + (t599 * t763 + (t676 * t683 - t691 * t735 - t683) * t605) * t795) * t772;
t578 = (t755 * t781 + (-t599 * t650 * t766 + t683 * (t599 * t792 - t741)) * t599) * t772;
t581 = (t590 * t795 - t781) * t804;
t596 = t599 ^ 2;
t602 = t605 ^ 2;
t656 = mrSges(3,1) * t691 + mrSges(3,2) * t697;
t632 = t681 * t656 * t692 - t732 * t683;
t729 = (-t572 * t769 + t578 * t632 + t581 * t673 + ((-mrSges(2,1) * t602 - t732 * (t602 + t596)) * t692 - 0.2e1 * t698 * t605 * (t656 * t599 + t605 * t797)) * t681 - t683 * t596 * t656) * t804;
t591 = t740 - t794;
t573 = (((t600 * t683 + t606 * t762) * t791 - (-pkin(5) * t606 + t743) * t734 + t683 * t591) * t606 - (-t600 * t762 + (-t678 * t683 + t693 * t734 + t683) * t606) * t794) * t771;
t579 = (t755 * t780 + (-t600 * t651 * t765 + t683 * (t600 * t791 - t740)) * t600) * t771;
t582 = (t591 * t794 - t780) * t803;
t597 = t600 ^ 2;
t603 = t606 ^ 2;
t657 = mrSges(3,1) * t693 + mrSges(3,2) * t699;
t633 = t681 * t657 * t694 - t731 * t683;
t728 = (-t573 * t768 + t579 * t633 + t582 * t673 + ((-mrSges(2,1) * t603 - t731 * (t603 + t597)) * t694 - 0.2e1 * t700 * t606 * (t657 * t600 + t606 * t797)) * t681 - t683 * t597 * t657) * t803;
t658 = -Ifges(3,5) * t689 - Ifges(3,6) * t695;
t684 = Ifges(3,1) - Ifges(3,2);
t727 = 0.2e1 * ((t598 * t799 + t684 * t776) * t695 + t779 * t798 + (t802 - 0.1e1) * t604 * Ifges(3,4)) * t598 - t580 * t770 + (t674 * t684 + t754 * t806 + t786) * t571 + t658 * t577;
t659 = -Ifges(3,5) * t691 - Ifges(3,6) * t697;
t726 = 0.2e1 * ((t599 * t799 + t684 * t775) * t697 + t778 * t798 + (t801 - 0.1e1) * t605 * Ifges(3,4)) * t599 - t581 * t769 + (t676 * t684 + t751 * t806 + t786) * t572 + t659 * t578;
t660 = -Ifges(3,5) * t693 - Ifges(3,6) * t699;
t725 = 0.2e1 * ((t600 * t799 + t684 * t774) * t699 + t777 * t798 + (t800 - 0.1e1) * t606 * Ifges(3,4)) * t600 - t582 * t768 + (t678 * t684 + t748 * t806 + t786) * t573 + t660 * t579;
t724 = Ifges(3,3) * t577 - t571 * t658 - t580 * t631 + t601 * (Ifges(3,4) * t802 + t684 * t754 - Ifges(3,4));
t723 = Ifges(3,3) * t578 - t572 * t659 - t581 * t632 + t602 * (Ifges(3,4) * t801 + t684 * t751 - Ifges(3,4));
t722 = Ifges(3,3) * t579 - t573 * t660 - t582 * t633 + t603 * (Ifges(3,4) * t800 + t684 * t748 - Ifges(3,4));
t718 = pkin(2) * t767 - t649 * t683;
t717 = pkin(2) * t766 - t650 * t683;
t716 = pkin(2) * t765 - t651 * t683;
t715 = t727 * t773;
t714 = t726 * t772;
t713 = t725 * t771;
t712 = t724 * t773;
t711 = t723 * t772;
t710 = t722 * t771;
t654 = pkin(5) * t694 + t700 * t787;
t653 = pkin(5) * t692 + t698 * t788;
t652 = pkin(5) * t690 + t696 * t789;
t615 = t654 * t682 + t716 * t680;
t614 = t653 * t682 + t717 * t680;
t613 = t652 * t682 + t718 * t680;
t1 = [-t725 * t630 * t737 - t726 * t629 * t738 - t727 * t628 * t739 + (t615 * t672 + t642 * t669) * t728 + (t614 * t671 + t641 * t668) * t729 + (t613 * t670 + t640 * t667) * t730 + (t724 * t619 * t739 + t723 * t620 * t738 + t722 * t621 * t737) * t706; t669 * t630 * t713 + t668 * t629 * t714 + t667 * t628 * t715 + (-t615 * t669 + t642 * t672) * t728 + (-t614 * t668 + t641 * t671) * t729 + (-t613 * t667 + t640 * t670) * t730 + (-t667 * t619 * t712 - t668 * t620 * t711 - t669 * t621 * t710) * t706; t627 * t713 + t626 * t714 + t625 * t715 + (-t654 * t680 + t716 * t682) * t728 + (-t653 * t680 + t717 * t682) * t729 + (-t652 * t680 + t718 * t682) * t730 + (-t622 * t712 - t623 * t711 - t624 * t710) * t706;];
taucX  = t1;
