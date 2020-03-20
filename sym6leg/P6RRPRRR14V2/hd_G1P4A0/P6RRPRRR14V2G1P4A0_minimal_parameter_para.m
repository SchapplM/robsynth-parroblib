% Return the minimum parameter vector for
% P6RRPRRR14V2G1P4A0
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d1,d2,d4,theta3]';
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
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% 
% Output:
% MPV [24x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-12 22:59
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = P6RRPRRR14V2G1P4A0_minimal_parameter_para(pkin, m, mrSges, Ifges, koppelP)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6),zeros(6,3)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6RRPRRR14V2G1P4A0_minimal_parameter_para: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6RRPRRR14V2G1P4A0_minimal_parameter_para: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P6RRPRRR14V2G1P4A0_minimal_parameter_para: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P6RRPRRR14V2G1P4A0_minimal_parameter_para: Ifges has to be [4x6] (double)'); 
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6RRPRRR14V2G1P4A0_minimal_parameter_para: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From minimal_parameter_parrob_matlab.m
unknown=NaN(24,1);
t1 = sin(pkin(4));
t2 = t1 ^ 2;
t3 = sin(pkin(10));
t4 = t3 ^ 2;
t5 = cos(pkin(5));
t6 = t5 ^ 2;
t8 = Ifges(3,1) * t6 * t4;
t9 = cos(pkin(10));
t10 = t3 * t9;
t13 = 0.2e1 * Ifges(3,4) * t6 * t10;
t14 = t5 * t3;
t15 = sin(pkin(5));
t18 = 0.2e1 * Ifges(3,5) * t15 * t14;
t19 = t9 ^ 2;
t21 = Ifges(3,2) * t6 * t19;
t22 = t5 * t9;
t25 = 0.2e1 * Ifges(3,6) * t15 * t22;
t26 = t15 ^ 2;
t27 = Ifges(3,3) * t26;
t33 = pkin(1) ^ 2;
t34 = pkin(8) ^ 2;
t49 = Ifges(3,1) * t19 - 0.2e1 * Ifges(3,4) * t10 + Ifges(3,2) * t4 + Ifges(2,1) - Ifges(2,2) - t13 + t18 - t21 + t25 - t27 - t8;
t52 = t19 - t4;
t62 = Ifges(3,1) * t15;
t67 = Ifges(3,2) * t15;
t73 = t15 * t5;
t77 = t6 - t26;
unknown(1,1) = Ifges(1,3) + (Ifges(2,2) + t8 + t13 - t18 + t21 - t25 + t27) * t2 + 0.2e1 * mrSges(2,3) * t2 * pkin(8) + m(2) * (t2 * t34 + t33);
unknown(2,1) = m(2) * pkin(1) + mrSges(1,1);
unknown(3,1) = -m(2) * t1 * pkin(8) - mrSges(2,3) * t1 + mrSges(1,2);
unknown(4,1) = t49;
unknown(5,1) = Ifges(3,1) * t5 * t10 + Ifges(3,4) * t5 * t52 - Ifges(3,5) * t15 * t9 - Ifges(3,2) * t5 * t10 + Ifges(3,6) * t15 * t3 + Ifges(2,4);
unknown(6,1) = Ifges(3,4) * t15 * t52 + Ifges(3,5) * t22 - Ifges(3,6) * t14 + t62 * t10 - t67 * t10 + Ifges(2,5);
unknown(7,1) = 0.2e1 * Ifges(3,4) * t73 * t10 + Ifges(3,5) * t77 * t3 + Ifges(3,6) * t77 * t9 + t67 * t5 * t19 + t62 * t5 * t4 - Ifges(3,3) * t73 + Ifges(2,6);
unknown(8,1) = Ifges(3,1) * t26 * t4 + 0.2e1 * Ifges(3,4) * t26 * t10 + Ifges(3,2) * t26 * t19 + Ifges(3,3) * t6 + Ifges(2,3) + t18 + t25;
unknown(9,1) = mrSges(2,1);
unknown(10,1) = mrSges(2,2);
unknown(11,1) = mrSges(3,1);
unknown(12,1) = mrSges(3,2);
unknown(13,1) = mrSges(3,3);
unknown(14,1) = m(3);
unknown(15,1) = Ifges(4,1);
unknown(16,1) = Ifges(4,4);
unknown(17,1) = Ifges(4,5);
unknown(18,1) = Ifges(4,2);
unknown(19,1) = Ifges(4,6);
unknown(20,1) = Ifges(4,3);
unknown(21,1) = mrSges(4,1);
unknown(22,1) = mrSges(4,2);
unknown(23,1) = mrSges(4,3);
unknown(24,1) = m(4);
MPV  = unknown;
