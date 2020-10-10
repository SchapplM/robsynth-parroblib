% Calculate Gravitation load for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_gravload_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:38
% EndTime: 2020-08-06 18:24:38
% DurationCPUTime: 0.39s
% Computational Cost: add. (423->104), mult. (567->179), div. (36->7), fcn. (375->18), ass. (0->77)
t658 = m(2) + m(3);
t657 = pkin(1) + pkin(5);
t651 = cos(qJ(3,3));
t687 = pkin(3) * t651;
t652 = cos(qJ(1,3));
t686 = pkin(3) * t652;
t653 = cos(qJ(3,2));
t685 = pkin(3) * t653;
t654 = cos(qJ(1,2));
t684 = pkin(3) * t654;
t655 = cos(qJ(3,1));
t683 = pkin(3) * t655;
t656 = cos(qJ(1,1));
t682 = pkin(3) * t656;
t681 = -mrSges(1,2) + mrSges(2,3);
t642 = legFrame(3,2);
t632 = sin(t642);
t635 = cos(t642);
t615 = t635 * g(1) - t632 * g(2);
t646 = sin(qJ(1,3));
t609 = g(3) * t646 - t652 * t615;
t612 = t632 * g(1) + t635 * g(2);
t645 = sin(qJ(3,3));
t639 = 0.1e1 / t645;
t680 = ((-t609 * mrSges(3,1) + mrSges(3,2) * t612) * t651 + t645 * (mrSges(3,1) * t612 + t609 * mrSges(3,2))) * t639;
t643 = legFrame(2,2);
t633 = sin(t643);
t636 = cos(t643);
t616 = t636 * g(1) - t633 * g(2);
t648 = sin(qJ(1,2));
t610 = g(3) * t648 - t654 * t616;
t613 = t633 * g(1) + t636 * g(2);
t647 = sin(qJ(3,2));
t640 = 0.1e1 / t647;
t679 = ((-t610 * mrSges(3,1) + mrSges(3,2) * t613) * t653 + t647 * (mrSges(3,1) * t613 + t610 * mrSges(3,2))) * t640;
t644 = legFrame(1,2);
t634 = sin(t644);
t637 = cos(t644);
t617 = t637 * g(1) - t634 * g(2);
t650 = sin(qJ(1,1));
t611 = g(3) * t650 - t656 * t617;
t614 = t634 * g(1) + t637 * g(2);
t649 = sin(qJ(3,1));
t641 = 0.1e1 / t649;
t678 = ((-t611 * mrSges(3,1) + mrSges(3,2) * t614) * t655 + t649 * (mrSges(3,1) * t614 + t611 * mrSges(3,2))) * t641;
t629 = t645 * pkin(3) + qJ(2,3);
t626 = 0.1e1 / t629;
t677 = t609 * t626;
t630 = t647 * pkin(3) + qJ(2,2);
t627 = 0.1e1 / t630;
t676 = t610 * t627;
t631 = t649 * pkin(3) + qJ(2,1);
t628 = 0.1e1 / t631;
t675 = t611 * t628;
t622 = m(2) * pkin(1) + t657 * m(3) + mrSges(1,1) - mrSges(2,2) + mrSges(3,3);
t621 = t622 * g(3);
t662 = t645 * mrSges(3,1) + mrSges(3,2) * t651 + t658 * qJ(2,3) + t681;
t674 = t626 * (t621 * t646 + (-t622 * t652 - t662 * t646) * t615 - t662 * t652 * g(3));
t661 = t647 * mrSges(3,1) + mrSges(3,2) * t653 + t658 * qJ(2,2) + t681;
t673 = t627 * (t621 * t648 + (-t622 * t654 - t661 * t648) * t616 - t661 * t654 * g(3));
t660 = t649 * mrSges(3,1) + mrSges(3,2) * t655 + t658 * qJ(2,1) + t681;
t672 = t628 * (t621 * t650 + (-t622 * t656 - t660 * t650) * t617 - t660 * t656 * g(3));
t671 = t651 * qJ(2,3);
t670 = t653 * qJ(2,2);
t669 = t655 * qJ(2,1);
t668 = t639 * t677;
t667 = t640 * t676;
t666 = t641 * t675;
t665 = t646 * t674;
t664 = t648 * t673;
t663 = t650 * t672;
t659 = 0.1e1 / pkin(3);
t638 = pkin(6) + t657;
t620 = qJ(2,1) * t656 - t638 * t650;
t619 = qJ(2,2) * t654 - t638 * t648;
t618 = qJ(2,3) * t652 - t638 * t646;
t1 = [t635 * t665 + t636 * t664 + t637 * t663 - g(1) * m(4) + (-t632 * t680 - t633 * t679 - t634 * t678) * t659 + (-((-t620 * t637 + t634 * t683) * t649 + (t655 - 0.1e1) * (t655 + 0.1e1) * t637 * t682 + t634 * t669) * t666 - ((-t619 * t636 + t633 * t685) * t647 + (t653 - 0.1e1) * (t653 + 0.1e1) * t636 * t684 + t633 * t670) * t667 - ((-t618 * t635 + t632 * t687) * t645 + (t651 - 0.1e1) * (t651 + 0.1e1) * t635 * t686 + t632 * t671) * t668) * t658; -t632 * t665 - t633 * t664 - t634 * t663 - g(2) * m(4) + (-t635 * t680 - t636 * t679 - t637 * t678) * t659 + (-((t620 * t634 + t637 * t683) * t649 + (-t655 ^ 2 + 0.1e1) * t634 * t682 + t637 * t669) * t666 - ((t619 * t633 + t636 * t685) * t647 + (-t653 ^ 2 + 0.1e1) * t633 * t684 + t636 * t670) * t667 - ((t618 * t632 + t635 * t687) * t645 + (-t651 ^ 2 + 0.1e1) * t632 * t686 + t635 * t671) * t668) * t658; t652 * t674 + t654 * t673 + t656 * t672 - g(3) * m(4) + (-(t650 * t631 + t638 * t656) * t675 - (t648 * t630 + t638 * t654) * t676 - (t646 * t629 + t638 * t652) * t677) * t658;];
taugX  = t1;
