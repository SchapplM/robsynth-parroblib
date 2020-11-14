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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:44
% EndTime: 2020-08-06 18:41:45
% DurationCPUTime: 1.76s
% Computational Cost: add. (3774->255), mult. (4833->430), div. (531->10), fcn. (3690->53), ass. (0->172)
t650 = sin(pkin(7));
t670 = -pkin(6) - pkin(5);
t612 = t670 * t650 - pkin(1);
t663 = cos(qJ(1,3));
t597 = t612 * t663;
t651 = cos(pkin(7));
t657 = sin(qJ(1,3));
t714 = t657 * t670;
t715 = t657 * t650;
t745 = t597 + pkin(2) * t715 - (pkin(2) * t663 - t714) * t651;
t665 = cos(qJ(1,2));
t598 = t612 * t665;
t659 = sin(qJ(1,2));
t710 = t659 * t670;
t711 = t659 * t650;
t744 = t598 + pkin(2) * t711 - (pkin(2) * t665 - t710) * t651;
t667 = cos(qJ(1,1));
t599 = t612 * t667;
t661 = sin(qJ(1,1));
t706 = t661 * t670;
t707 = t661 * t650;
t743 = t599 + pkin(2) * t707 - (pkin(2) * t667 - t706) * t651;
t629 = t651 * pkin(1);
t620 = t629 + pkin(2);
t742 = 0.2e1 * pkin(1);
t741 = 0.2e1 * pkin(2);
t740 = 0.2e1 * t670;
t739 = pkin(1) * t650;
t662 = cos(qJ(3,3));
t621 = t662 * pkin(3) + pkin(2);
t664 = cos(qJ(3,2));
t622 = t664 * pkin(3) + pkin(2);
t666 = cos(qJ(3,1));
t623 = t666 * pkin(3) + pkin(2);
t656 = sin(qJ(3,3));
t738 = t656 * mrSges(3,2);
t658 = sin(qJ(3,2));
t737 = t658 * mrSges(3,2);
t660 = sin(qJ(3,1));
t736 = t660 * mrSges(3,2);
t672 = 0.1e1 / pkin(3);
t720 = t621 * t663;
t735 = ((-t714 + t720) * t651 - t597 - t621 * t715) * t672;
t719 = t622 * t665;
t734 = ((-t710 + t719) * t651 - t598 - t622 * t711) * t672;
t718 = t623 * t667;
t733 = ((-t706 + t718) * t651 - t599 - t623 * t707) * t672;
t640 = qJ(1,3) + pkin(7);
t626 = cos(t640);
t697 = pkin(7) + qJ(3,3);
t700 = -pkin(7) + qJ(3,3);
t732 = (t626 * t740 + sin(t640) * t741 + t657 * t742 + (sin(qJ(1,3) - t700) + sin(qJ(1,3) + t697)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,3)) + t656 * t741 + (sin(t697) + sin(t700)) * pkin(1));
t641 = qJ(1,2) + pkin(7);
t627 = cos(t641);
t698 = pkin(7) + qJ(3,2);
t701 = -pkin(7) + qJ(3,2);
t731 = (t627 * t740 + sin(t641) * t741 + t659 * t742 + (sin(qJ(1,2) - t701) + sin(qJ(1,2) + t698)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,2)) + t658 * t741 + (sin(t698) + sin(t701)) * pkin(1));
t642 = qJ(1,1) + pkin(7);
t628 = cos(t642);
t699 = pkin(7) + qJ(3,1);
t702 = -pkin(7) + qJ(3,1);
t730 = (t628 * t740 + sin(t642) * t741 + t661 * t742 + (sin(qJ(1,1) - t702) + sin(qJ(1,1) + t699)) * pkin(3)) / (pkin(3) * sin(0.2e1 * qJ(3,1)) + t660 * t741 + (sin(t699) + sin(t702)) * pkin(1));
t653 = legFrame(3,2);
t614 = t653 + t640;
t615 = -t653 + t640;
t729 = -sin(t614) / 0.2e1 - sin(t615) / 0.2e1;
t654 = legFrame(2,2);
t616 = t654 + t641;
t617 = -t654 + t641;
t728 = -sin(t616) / 0.2e1 - sin(t617) / 0.2e1;
t655 = legFrame(1,2);
t618 = t655 + t642;
t619 = -t655 + t642;
t727 = -sin(t618) / 0.2e1 - sin(t619) / 0.2e1;
t726 = cos(t615) / 0.2e1 - cos(t614) / 0.2e1;
t725 = cos(t617) / 0.2e1 - cos(t616) / 0.2e1;
t724 = cos(t619) / 0.2e1 - cos(t618) / 0.2e1;
t600 = 0.1e1 / (t629 + t621);
t644 = 0.1e1 / t656;
t723 = t600 * t644;
t601 = 0.1e1 / (t629 + t622);
t645 = 0.1e1 / t658;
t722 = t601 * t645;
t602 = 0.1e1 / (t629 + t623);
t646 = 0.1e1 / t660;
t721 = t602 * t646;
t630 = sin(t653);
t717 = t656 * t630;
t633 = cos(t653);
t716 = t656 * t633;
t631 = sin(t654);
t713 = t658 * t631;
t634 = cos(t654);
t712 = t658 * t634;
t632 = sin(t655);
t709 = t660 * t632;
t635 = cos(t655);
t708 = t660 * t635;
t705 = t620 * mrSges(3,1);
t704 = -0.2e1 * mrSges(3,2) * pkin(2);
t703 = -0.2e1 * t629;
t647 = t662 ^ 2;
t696 = pkin(3) * (-t663 * t651 + t715) * t647;
t648 = t664 ^ 2;
t695 = pkin(3) * (-t665 * t651 + t711) * t648;
t649 = t666 ^ 2;
t694 = pkin(3) * (-t667 * t651 + t707) * t649;
t693 = -m(3) * pkin(2) - mrSges(2,1);
t692 = ((t621 * t657 + t663 * t670) * t651 - t612 * t657 + t650 * t720) * t644 * t662;
t691 = t630 * t735;
t690 = t633 * t735;
t689 = ((t622 * t659 + t665 * t670) * t651 - t612 * t659 + t650 * t719) * t645 * t664;
t688 = t631 * t734;
t687 = t634 * t734;
t686 = ((t623 * t661 + t667 * t670) * t651 - t612 * t661 + t650 * t718) * t646 * t666;
t685 = t632 * t733;
t684 = t635 * t733;
t683 = t672 * t732;
t682 = t672 * t731;
t681 = t672 * t730;
t680 = pkin(5) + t739;
t607 = t680 * mrSges(3,1) - Ifges(3,5);
t608 = -t680 * mrSges(3,2) + Ifges(3,6);
t582 = -t656 * t607 + t608 * t662;
t679 = t582 * t644 * t735;
t583 = -t658 * t607 + t608 * t664;
t678 = t583 * t645 * t734;
t584 = -t660 * t607 + t608 * t666;
t677 = t584 * t646 * t733;
t671 = m(2) + m(3);
t673 = Ifges(3,1) + Ifges(1,3) + Ifges(2,3) + 0.2e1 * (m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3)) * t739 + (pkin(2) ^ 2 + pkin(5) ^ 2) * m(3) + t671 * pkin(1) ^ 2 + 0.2e1 * mrSges(3,3) * pkin(5);
t652 = -Ifges(3,1) + Ifges(3,2);
t611 = t666 * mrSges(3,1) - t736;
t610 = t664 * mrSges(3,1) - t737;
t609 = t662 * mrSges(3,1) - t738;
t572 = t652 * t649 + 0.2e1 * (Ifges(3,4) * t660 + t705) * t666 + (t693 + t736) * t703 + t660 * t704 + t673;
t571 = t652 * t648 + 0.2e1 * (Ifges(3,4) * t658 + t705) * t664 + (t693 + t737) * t703 + t658 * t704 + t673;
t570 = t652 * t647 + 0.2e1 * (Ifges(3,4) * t656 + t705) * t662 + (t693 + t738) * t703 + t656 * t704 + t673;
t569 = -t635 * t694 + (pkin(3) * t709 - t743 * t635) * t666 + t620 * t709;
t568 = -t634 * t695 + (pkin(3) * t713 - t744 * t634) * t664 + t620 * t713;
t567 = -t633 * t696 + (pkin(3) * t717 - t745 * t633) * t662 + t620 * t717;
t566 = t632 * t694 + (pkin(3) * t708 + t743 * t632) * t666 + t620 * t708;
t565 = t631 * t695 + (pkin(3) * t712 + t744 * t631) * t664 + t620 * t712;
t564 = t630 * t696 + (pkin(3) * t716 + t745 * t630) * t662 + t620 * t716;
t563 = -t602 * t671 * t686 + t611 * t681;
t562 = -t601 * t671 * t689 + t610 * t682;
t561 = -t600 * t671 * t692 + t609 * t683;
t560 = -t628 * t602 * t572 + t584 * t681;
t559 = -t627 * t601 * t571 + t583 * t682;
t558 = -t626 * t600 * t570 + t582 * t683;
t557 = (t569 * t671 - t611 * t684) * t721;
t556 = (t568 * t671 - t610 * t687) * t722;
t555 = (t567 * t671 - t609 * t690) * t723;
t554 = (t566 * t671 + t611 * t685) * t721;
t553 = (t565 * t671 + t610 * t688) * t722;
t552 = (t564 * t671 + t609 * t691) * t723;
t551 = Ifges(3,3) * t681 + (-t584 * t628 - t611 * t686) * t602;
t550 = Ifges(3,3) * t682 + (-t583 * t627 - t610 * t689) * t601;
t549 = Ifges(3,3) * t683 + (-t582 * t626 - t609 * t692) * t600;
t548 = (t572 * t724 + t632 * t677) * t602;
t547 = (t571 * t725 + t631 * t678) * t601;
t546 = (t570 * t726 + t630 * t679) * t600;
t545 = (t572 * t727 - t635 * t677) * t602;
t544 = (t571 * t728 - t634 * t678) * t601;
t543 = (t570 * t729 - t633 * t679) * t600;
t542 = (t584 * t724 + (Ifges(3,3) * t685 + t566 * t611) * t646) * t602;
t541 = (t583 * t725 + (Ifges(3,3) * t688 + t565 * t610) * t645) * t601;
t540 = (t582 * t726 + (Ifges(3,3) * t691 + t564 * t609) * t644) * t600;
t539 = (t584 * t727 + (-Ifges(3,3) * t684 + t569 * t611) * t646) * t602;
t538 = (t583 * t728 + (-Ifges(3,3) * t687 + t568 * t610) * t645) * t601;
t537 = (t582 * t729 + (-Ifges(3,3) * t690 + t567 * t609) * t644) * t600;
t1 = [m(4) + (t545 * t727 + (-t539 * t684 + t557 * t569) * t646) * t602 + (t544 * t728 + (-t538 * t687 + t556 * t568) * t645) * t601 + (t543 * t729 + (-t537 * t690 + t555 * t567) * t644) * t600, (t545 * t724 + (t539 * t685 + t557 * t566) * t646) * t602 + (t544 * t725 + (t538 * t688 + t556 * t565) * t645) * t601 + (t543 * t726 + (t537 * t691 + t555 * t564) * t644) * t600, (-t545 * t628 - t557 * t686) * t602 + (-t544 * t627 - t556 * t689) * t601 + (-t543 * t626 - t555 * t692) * t600 + (t537 * t732 + t538 * t731 + t539 * t730) * t672; (t548 * t727 + (-t542 * t684 + t554 * t569) * t646) * t602 + (t547 * t728 + (-t541 * t687 + t553 * t568) * t645) * t601 + (t546 * t729 + (-t540 * t690 + t552 * t567) * t644) * t600, m(4) + (t548 * t724 + (t542 * t685 + t554 * t566) * t646) * t602 + (t547 * t725 + (t541 * t688 + t553 * t565) * t645) * t601 + (t546 * t726 + (t540 * t691 + t552 * t564) * t644) * t600, (-t548 * t628 - t554 * t686) * t602 + (-t547 * t627 - t553 * t689) * t601 + (-t546 * t626 - t552 * t692) * t600 + (t540 * t732 + t541 * t731 + t542 * t730) * t672; (t560 * t727 + (-t551 * t684 + t563 * t569) * t646) * t602 + (t559 * t728 + (-t550 * t687 + t562 * t568) * t645) * t601 + (t558 * t729 + (-t549 * t690 + t561 * t567) * t644) * t600, (t560 * t724 + (t551 * t685 + t563 * t566) * t646) * t602 + (t559 * t725 + (t550 * t688 + t562 * t565) * t645) * t601 + (t558 * t726 + (t549 * t691 + t561 * t564) * t644) * t600, m(4) + (-t560 * t628 - t563 * t686) * t602 + (-t559 * t627 - t562 * t689) * t601 + (-t558 * t626 - t561 * t692) * t600 + (t549 * t732 + t550 * t731 + t551 * t730) * t672;];
MX  = t1;
