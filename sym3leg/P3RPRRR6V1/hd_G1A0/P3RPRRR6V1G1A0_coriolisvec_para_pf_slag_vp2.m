% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPRRR6V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:04
% EndTime: 2020-08-06 18:32:06
% DurationCPUTime: 2.38s
% Computational Cost: add. (11718->270), mult. (10161->369), div. (1347->13), fcn. (5877->68), ass. (0->179)
t776 = sin(pkin(7)) * pkin(1);
t804 = pkin(5) + t776;
t809 = t804 * mrSges(3,2);
t678 = cos(pkin(7)) * pkin(1);
t670 = t678 + pkin(2);
t808 = mrSges(3,1) * t670;
t701 = Ifges(3,1) - Ifges(3,2);
t704 = sin(qJ(3,1));
t707 = cos(qJ(3,1));
t694 = t707 ^ 2;
t785 = -0.2e1 * t694;
t803 = mrSges(3,2) * t670;
t807 = (t785 + 0.1e1) * Ifges(3,4) + (-t701 * t704 + t803) * t707;
t703 = sin(qJ(3,2));
t706 = cos(qJ(3,2));
t693 = t706 ^ 2;
t786 = -0.2e1 * t693;
t806 = (t786 + 0.1e1) * Ifges(3,4) + (-t701 * t703 + t803) * t706;
t702 = sin(qJ(3,3));
t705 = cos(qJ(3,3));
t692 = t705 ^ 2;
t787 = -0.2e1 * t692;
t805 = (t787 + 0.1e1) * Ifges(3,4) + (-t701 * t702 + t803) * t705;
t754 = 0.2e1 * pkin(2);
t672 = qJ(1,3) + legFrame(3,3);
t661 = pkin(7) + t672;
t640 = sin(t661);
t643 = cos(t661);
t648 = qJ(3,3) + t661;
t649 = -qJ(3,3) + t661;
t762 = sin(t648) + sin(t649);
t713 = -pkin(6) - pkin(5);
t778 = -0.2e1 * t713;
t788 = -0.2e1 * pkin(2);
t789 = -0.2e1 * pkin(1);
t588 = t643 * t778 + sin(t672) * t789 + t640 * t788 - t762 * pkin(3);
t759 = cos(t648) + cos(t649);
t777 = 0.2e1 * t713;
t591 = t640 * t777 + cos(t672) * t789 + t643 * t788 - t759 * pkin(3);
t716 = 0.2e1 * qJ(3,3);
t688 = sin(t716);
t675 = pkin(3) * t688;
t784 = 0.2e1 * t702;
t751 = pkin(2) * t784;
t606 = 0.1e1 / (t751 + t675 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1));
t711 = xDP(2);
t712 = xDP(1);
t720 = 0.1e1 / pkin(3);
t579 = (t588 * t711 + t591 * t712) * t720 * t606;
t802 = 0.2e1 * t579;
t673 = qJ(1,2) + legFrame(2,3);
t662 = pkin(7) + t673;
t641 = sin(t662);
t644 = cos(t662);
t652 = qJ(3,2) + t662;
t653 = -qJ(3,2) + t662;
t761 = sin(t652) + sin(t653);
t589 = t644 * t778 + sin(t673) * t789 + t641 * t788 - t761 * pkin(3);
t758 = cos(t652) + cos(t653);
t592 = t641 * t777 + cos(t673) * t789 + t644 * t788 - t758 * pkin(3);
t717 = 0.2e1 * qJ(3,2);
t689 = sin(t717);
t676 = pkin(3) * t689;
t783 = 0.2e1 * t703;
t750 = pkin(2) * t783;
t607 = 0.1e1 / (t750 + t676 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1));
t580 = (t589 * t711 + t592 * t712) * t720 * t607;
t801 = 0.2e1 * t580;
t674 = qJ(1,1) + legFrame(1,3);
t663 = pkin(7) + t674;
t642 = sin(t663);
t645 = cos(t663);
t656 = qJ(3,1) + t663;
t657 = -qJ(3,1) + t663;
t760 = sin(t656) + sin(t657);
t590 = t645 * t778 + sin(t674) * t789 + t642 * t788 - t760 * pkin(3);
t757 = cos(t656) + cos(t657);
t593 = t642 * t777 + cos(t674) * t789 + t645 * t788 - t757 * pkin(3);
t718 = 0.2e1 * qJ(3,1);
t690 = sin(t718);
t677 = pkin(3) * t690;
t782 = 0.2e1 * t704;
t749 = pkin(2) * t782;
t608 = 0.1e1 / (t749 + t677 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1));
t581 = (t590 * t711 + t593 * t712) * t720 * t608;
t800 = 0.2e1 * t581;
t578 = t581 ^ 2;
t773 = t707 * pkin(3) + pkin(2);
t619 = t678 + t773;
t615 = 0.1e1 / t619;
t599 = (-t642 * t712 + t645 * t711) * t615;
t596 = t599 ^ 2;
t719 = pkin(3) ^ 2;
t721 = pkin(2) ^ 2;
t722 = pkin(1) ^ 2;
t753 = -0.2e1 * t776;
t734 = t721 + t722 + (t753 + t713) * t713;
t752 = 0.2e1 * t678;
t609 = pkin(2) * t752 + t719 / 0.2e1 + t734;
t764 = t704 * t670;
t612 = 0.1e1 / (t677 + 0.2e1 * t764);
t660 = pkin(6) + t804;
t779 = -0.2e1 * t707;
t566 = (t596 * t609 * t779 + (-0.2e1 * t578 * t619 + ((-cos(t718) * t619 + t670 * t785 - t670) * t599 + (t707 * t782 + t690) * t660 * t581) * t599) * pkin(3)) * t612;
t745 = t581 * t704 * pkin(3);
t763 = pkin(3) * t754;
t572 = (-t660 * t745 + (t694 * t719 + t707 * t763 + t752 * t773 + t734) * t599) / (t677 / 0.2e1 + t764) * t599 * t720 + (t707 * t581 * t619 - t599 * t660 * t704) * t612 * t800;
t683 = t704 * mrSges(3,2);
t622 = -t707 * mrSges(3,1) + t683;
t715 = m(2) + m(3);
t790 = -t715 * t566 + t622 * t572 - t578 * (mrSges(3,1) * t704 + mrSges(3,2) * t707);
t798 = t790 * t608 / 0.2e1;
t577 = t580 ^ 2;
t774 = t706 * pkin(3) + pkin(2);
t618 = t678 + t774;
t614 = 0.1e1 / t618;
t598 = (-t641 * t712 + t644 * t711) * t614;
t595 = t598 ^ 2;
t765 = t703 * t670;
t611 = 0.1e1 / (t676 + 0.2e1 * t765);
t780 = -0.2e1 * t706;
t565 = (t595 * t609 * t780 + (-0.2e1 * t577 * t618 + ((-cos(t717) * t618 + t670 * t786 - t670) * t598 + (t706 * t783 + t689) * t660 * t580) * t598) * pkin(3)) * t611;
t746 = t580 * t703 * pkin(3);
t571 = (-t660 * t746 + (t693 * t719 + t706 * t763 + t752 * t774 + t734) * t598) / (t676 / 0.2e1 + t765) * t598 * t720 + (t706 * t580 * t618 - t598 * t660 * t703) * t611 * t801;
t682 = t703 * mrSges(3,2);
t624 = -t706 * mrSges(3,1) + t682;
t791 = -t715 * t565 + t624 * t571 - t577 * (mrSges(3,1) * t703 + mrSges(3,2) * t706);
t797 = t791 * t607 / 0.2e1;
t576 = t579 ^ 2;
t775 = t705 * pkin(3) + pkin(2);
t617 = t678 + t775;
t613 = 0.1e1 / t617;
t597 = (-t640 * t712 + t643 * t711) * t613;
t594 = t597 ^ 2;
t766 = t702 * t670;
t610 = 0.1e1 / (t675 + 0.2e1 * t766);
t781 = -0.2e1 * t705;
t564 = (t594 * t609 * t781 + (-0.2e1 * t576 * t617 + ((-cos(t716) * t617 + t670 * t787 - t670) * t597 + (t705 * t784 + t688) * t660 * t579) * t597) * pkin(3)) * t610;
t747 = t579 * t702 * pkin(3);
t570 = (-t660 * t747 + (t692 * t719 + t705 * t763 + t752 * t775 + t734) * t597) / (t675 / 0.2e1 + t766) * t597 * t720 + (t705 * t579 * t617 - t597 * t660 * t702) * t610 * t802;
t681 = t702 * mrSges(3,2);
t623 = -t705 * mrSges(3,1) + t681;
t792 = -t715 * t564 + t623 * t570 - t576 * (mrSges(3,1) * t702 + mrSges(3,2) * t705);
t796 = t792 * t606 / 0.2e1;
t739 = -t713 + t776;
t574 = (t597 * t739 - 0.2e1 * t747) * t613 * t597;
t680 = mrSges(3,1) * pkin(5) - Ifges(3,5);
t748 = mrSges(3,1) * t776;
t620 = t748 + t680;
t621 = Ifges(3,6) - t809;
t603 = t702 * t620 - t621 * t705;
t723 = -(pkin(5) ^ 2 + t721) * m(3) - t715 * t722 - Ifges(3,1) - Ifges(1,3) - Ifges(2,3) + (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t753 - 0.2e1 * mrSges(3,3) * pkin(5);
t725 = t748 / 0.2e1 + t680 / 0.2e1;
t726 = Ifges(3,6) / 0.2e1 - t809 / 0.2e1;
t744 = -m(3) * pkin(2) - mrSges(2,1);
t795 = t613 * (((t726 * t702 + t725 * t705) * t579 + (t702 * t808 + t805) * t597) * t802 - (t701 * t692 + (Ifges(3,4) * t702 + t808) * t781 + (t681 + t744) * t752 + mrSges(3,2) * t751 + t723) * t574 - t603 * t570);
t575 = (t598 * t739 - 0.2e1 * t746) * t614 * t598;
t604 = t703 * t620 - t621 * t706;
t794 = t614 * (((t726 * t703 + t725 * t706) * t580 + (t703 * t808 + t806) * t598) * t801 - (t701 * t693 + (Ifges(3,4) * t703 + t808) * t780 + (t682 + t744) * t752 + mrSges(3,2) * t750 + t723) * t575 - t604 * t571);
t573 = (t599 * t739 - 0.2e1 * t745) * t615 * t599;
t605 = t704 * t620 - t621 * t707;
t793 = t615 * (((t726 * t704 + t725 * t707) * t581 + (t704 * t808 + t807) * t599) * t800 - (t701 * t694 + (Ifges(3,4) * t704 + t808) * t779 + (t683 + t744) * t752 + mrSges(3,2) * t749 + t723) * t573 - t605 * t572);
t755 = 0.2e1 * pkin(1);
t732 = (-Ifges(3,3) * t570 + t623 * t564 + t603 * t574 + t594 * (mrSges(3,1) * t766 + t805)) * t606;
t731 = (-Ifges(3,3) * t571 + t624 * t565 + t604 * t575 + t595 * (mrSges(3,1) * t765 + t806)) * t607;
t730 = (-Ifges(3,3) * t572 + t622 * t566 + t605 * t573 + t596 * (mrSges(3,1) * t764 + t807)) * t608;
t669 = -qJ(3,1) + t674;
t668 = qJ(3,1) + t674;
t667 = -qJ(3,2) + t673;
t666 = qJ(3,2) + t673;
t665 = -qJ(3,3) + t672;
t664 = qJ(3,3) + t672;
t658 = -0.2e1 * qJ(3,1) + t663;
t655 = t718 + t663;
t654 = -0.2e1 * qJ(3,2) + t662;
t651 = t717 + t662;
t650 = -0.2e1 * qJ(3,3) + t661;
t647 = t716 + t661;
t1 = [t642 * t793 + t641 * t794 + t640 * t795 + (t591 * t732 + t592 * t731 + t593 * t730) * t720 + (t757 * t754 + (cos(t669) + cos(t668)) * t755 + t760 * t778 + (cos(t658) + cos(t655) + 0.2e1 * t645) * pkin(3)) * t798 + (t758 * t754 + (cos(t667) + cos(t666)) * t755 + t761 * t778 + (cos(t654) + cos(t651) + 0.2e1 * t644) * pkin(3)) * t797 + (t759 * t754 + (cos(t665) + cos(t664)) * t755 + t762 * t778 + (cos(t650) + cos(t647) + 0.2e1 * t643) * pkin(3)) * t796; -t645 * t793 - t644 * t794 - t643 * t795 + (t588 * t732 + t589 * t731 + t590 * t730) * t720 + (t760 * t754 + (sin(t669) + sin(t668)) * t755 + t757 * t777 + (sin(t658) + sin(t655) + 0.2e1 * t642) * pkin(3)) * t798 + (t761 * t754 + (sin(t667) + sin(t666)) * t755 + t758 * t777 + (sin(t654) + sin(t651) + 0.2e1 * t641) * pkin(3)) * t797 + (t762 * t754 + (sin(t665) + sin(t664)) * t755 + t759 * t777 + (sin(t650) + sin(t647) + 0.2e1 * t640) * pkin(3)) * t796; t790 + t791 + t792;];
taucX  = t1;
