% Calculate Gravitation load for parallel robot
% P3RPRR1G3P3A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:27
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRR1G3P3A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G3P3A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:26:53
% EndTime: 2020-03-09 21:26:53
% DurationCPUTime: 0.30s
% Computational Cost: add. (579->132), mult. (564->133), div. (27->4), fcn. (330->66), ass. (0->87)
t663 = m(2) + m(3);
t693 = mrSges(3,1) * g(3);
t692 = mrSges(3,2) * g(3);
t691 = g(3) * mrSges(1,2);
t690 = g(3) * mrSges(2,2);
t656 = legFrame(3,2);
t646 = sin(t656);
t649 = cos(t656);
t595 = t649 * g(1) - t646 * g(2);
t671 = pkin(7) + qJ(3,3);
t637 = qJ(1,3) + t671;
t627 = sin(t637);
t583 = (mrSges(3,2) * t595 + t693) * cos(t637) + t627 * (mrSges(3,1) * t595 - t692);
t598 = 0.1e1 / (pkin(1) * sin(t671) + sin(qJ(3,3)) * pkin(2));
t689 = t583 * t598;
t657 = legFrame(2,2);
t647 = sin(t657);
t650 = cos(t657);
t596 = t650 * g(1) - t647 * g(2);
t672 = pkin(7) + qJ(3,2);
t638 = qJ(1,2) + t672;
t628 = sin(t638);
t584 = (mrSges(3,2) * t596 + t693) * cos(t638) + t628 * (mrSges(3,1) * t596 - t692);
t599 = 0.1e1 / (pkin(1) * sin(t672) + sin(qJ(3,2)) * pkin(2));
t688 = t584 * t599;
t658 = legFrame(1,2);
t648 = sin(t658);
t651 = cos(t658);
t597 = t651 * g(1) - t648 * g(2);
t673 = pkin(7) + qJ(3,1);
t639 = qJ(1,1) + t673;
t629 = sin(t639);
t585 = (mrSges(3,2) * t597 + t693) * cos(t639) + t629 * (mrSges(3,1) * t597 - t692);
t600 = 0.1e1 / (pkin(1) * sin(t673) + sin(qJ(3,1)) * pkin(2));
t687 = t585 * t600;
t686 = (t646 * g(1) + t649 * g(2)) * t663;
t685 = (t647 * g(1) + t650 * g(2)) * t663;
t684 = (t648 * g(1) + t651 * g(2)) * t663;
t620 = t663 * pkin(1) + mrSges(1,1);
t613 = g(3) * t620;
t652 = m(3) * pkin(2) + mrSges(2,1);
t633 = g(3) * t652;
t653 = qJ(1,3) + pkin(7);
t634 = sin(t653);
t659 = sin(qJ(1,3));
t683 = t598 * ((t595 * mrSges(2,2) + t633) * cos(t653) + (t652 * t595 - t690) * t634 + (mrSges(1,2) * t595 + t613) * cos(qJ(1,3)) + t659 * (t620 * t595 - t691) + t583);
t654 = qJ(1,2) + pkin(7);
t635 = sin(t654);
t660 = sin(qJ(1,2));
t682 = t599 * ((t596 * mrSges(2,2) + t633) * cos(t654) + (t652 * t596 - t690) * t635 + (mrSges(1,2) * t596 + t613) * cos(qJ(1,2)) + t660 * (t620 * t596 - t691) + t584);
t655 = qJ(1,1) + pkin(7);
t636 = sin(t655);
t661 = sin(qJ(1,1));
t681 = t600 * ((t597 * mrSges(2,2) + t633) * cos(t655) + (t652 * t597 - t690) * t636 + (mrSges(1,2) * t597 + t613) * cos(qJ(1,1)) + t661 * (t620 * t597 - t691) + t585);
t664 = 0.1e1 / pkin(3);
t680 = t664 / 0.2e1;
t621 = t656 + t653;
t614 = qJ(3,3) + t621;
t622 = -t656 + t653;
t615 = qJ(3,3) + t622;
t679 = -sin(t614) + sin(t615);
t623 = t657 + t654;
t616 = qJ(3,2) + t623;
t624 = -t657 + t654;
t617 = qJ(3,2) + t624;
t678 = -sin(t616) + sin(t617);
t625 = t658 + t655;
t618 = qJ(3,1) + t625;
t626 = -t658 + t655;
t619 = qJ(3,1) + t626;
t677 = -sin(t618) + sin(t619);
t676 = cos(t615) + cos(t614);
t675 = cos(t617) + cos(t616);
t674 = cos(t619) + cos(t618);
t670 = t683 / 0.2e1;
t669 = t682 / 0.2e1;
t668 = t681 / 0.2e1;
t667 = t680 * t689;
t666 = t680 * t688;
t665 = t680 * t687;
t645 = qJ(1,1) - t658;
t644 = qJ(1,1) + t658;
t643 = qJ(1,2) - t657;
t642 = qJ(1,2) + t657;
t641 = qJ(1,3) - t656;
t640 = qJ(1,3) + t656;
t1 = [t674 * t668 - t648 * t684 + (-t674 * pkin(3) + (-cos(t626) - cos(t625)) * pkin(2) + (-cos(t645) - cos(t644)) * pkin(1)) * t665 + t675 * t669 - t647 * t685 + (-t675 * pkin(3) + (-cos(t624) - cos(t623)) * pkin(2) + (-cos(t643) - cos(t642)) * pkin(1)) * t666 + t676 * t670 - t646 * t686 + (-t676 * pkin(3) + (-cos(t622) - cos(t621)) * pkin(2) + (-cos(t641) - cos(t640)) * pkin(1)) * t667 - g(1) * m(4); t677 * t668 - t651 * t684 - (t677 * pkin(3) + (-sin(t625) + sin(t626)) * pkin(2) + (-sin(t644) + sin(t645)) * pkin(1)) * t665 + t678 * t669 - t650 * t685 - (t678 * pkin(3) + (-sin(t623) + sin(t624)) * pkin(2) + (-sin(t642) + sin(t643)) * pkin(1)) * t666 + t679 * t670 - t649 * t686 - (t679 * pkin(3) + (-sin(t621) + sin(t622)) * pkin(2) + (-sin(t640) + sin(t641)) * pkin(1)) * t667 - g(2) * m(4); -t627 * t683 - t628 * t682 - t629 * t681 - g(3) * m(4) + ((t661 * pkin(1) + pkin(2) * t636 + pkin(3) * t629) * t687 + (t660 * pkin(1) + pkin(2) * t635 + pkin(3) * t628) * t688 + (t659 * pkin(1) + pkin(2) * t634 + pkin(3) * t627) * t689) * t664;];
taugX  = t1;
